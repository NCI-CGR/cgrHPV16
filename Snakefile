# vim: ft=python
import os
import collections
import glob
import itertools
import pandas
import pysam
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq

configfile: 'hpv_config.yaml'
workdir: os.environ['PWD']
shell.executable('bash')

hpv_ref = config['hpv_reference'] # currently, all bams are initially mapped to alpha reference
hpv_types = sorted(config['hpvtypes']) # just the numbers


# NOTE - OLD BUILDS of HPV16 may have a chomosome name that is uncompatible with this pipeline.

cov_dev = config['cov_dev'] # this is used for glu coverage script
GENES = config['genes'] # list of genes in gtf

bed = {} # e.g. bed['HPV16_Ref'] = 7906
for hpv in hpv_types:
    len_file = open(config['len_bed'] %hpv, 'r')
    for line in len_file:
        (type, length) = (line.split()[0], int(line.split()[2]))
        if config['padding'] == True:
            bed[type] = length - 400
        else:
            bed[type] = length
    len_file.close()

## ---- The parser will have to be customized for each run ---- ##
def parse_sampleID(filename):
    base = filename.split('/')[-1] 
    # For filenames IonXpress_046_SC075985.tmap.bam
    if base.startswith('I'):
        return base.rsplit('_', 3)[0].split('_', 2)[2]
    # otherwise PAP1111_2001_IonXpress_20.tmap.bam
    else:
        return base.split('_')[0]

bamfiles = sorted(glob.glob(config['tmap_path']), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(bamfiles, parse_sampleID):
    d[key] = list(value)

# We now have a dictionary with the sample ID as the keys and a list of 
# paths to the bam(s) as the value.
# e.g. d = {sample01: ['a/sample01.bam', 'b/sample01.bam'],
#           sample02: ['a/sample02.bam', 'c/sample02.bam'],
#           sample03: ['d/sample03.bam']}

sampleIDs = d.keys()

TARGETS =   [expand('bams/{sampleID}.bam', sampleID=sampleIDs),
            expand('reports/fasta/%s.HPV{hpv_type}.N-%d.fasta' %(config['deliver_proj'], config['fasta_n']), hpv_type=hpv_types),
            'reports/type_summary.tsv'
            ]

# eventually work these into the config yaml
include: 'qc_Snakefile' # creates the fastqc summaries and multiqc report
TARGETS += ['reports/filtered_read_count.tsv']
include: 'annotation_Snakefile' # annotates vcf and creates snpeff multiqc report
TARGETS += ['multiqc/snpeff_report.html', 'reports/%s_all_vcf_tables.txt' %config['deliver_proj']]
TARGETS += ['reports/median_gene_coverage.txt']


# These rules run on the host node and are not submitted to the cluster.
localrules: all

#--------------------------------------------------------------------------
rule all:
    input: TARGETS

def link_rule_input_files(wildcards):
    return d[wildcards.sampleID]

#--------------------------------------------------------------------------
rule link:
    input: link_rule_input_files
    output: 'bams/{sampleID}.bam'
    params:
        bam = 'temp/{sampleID}/{sampleID}.temp.bam',
        sam = 'temp/{sampleID}/{sampleID}.temp.sam'
    run:
        if (len(input) > 1):
            print(os.path.dirname(params.bam))
            shell('mkdir -p %s' %os.path.dirname(params.bam))
            shell('samtools merge {params.bam} {input}')

            # all SM fields in @RG must be identical
            # create a samfile with the fixed header
            sam = pysam.Samfile(params.bam, 'rb')
            header = sam.header
            for i in header['RG']:
                i['SM'] = wildcards.sampleID
            outfile = pysam.Samfile(params.sam, 'wh', header=header)

            # add the reads from the original merged bam
            shell('samtools view {params.bam} >> {params.sam}')

            # convert back to bam
            shell('samtools view -h -b {params.sam} > {output}')
            shell('rm {params.bam} {params.sam}') 
        else:
            shell('cd bams; ln -s ../{input} {wildcards.sampleID}.bam && touch -h {wildcards.sampleID}.bam')

#--------------------------------------------------------------------------
rule mapq_filter:
    input: rules.link.output
    output: 'mapq_filter/{sampleID}.filtered.bam'
    threads: 2
    params: 
        mapq = int(config["aq_filter"]),
        temp = 'temp/{sampleID}/{sampleID}.aq.temp',
        pre =  'temp/{sampleID}/{sampleID}.sort.temp'
    run:
        shell('mkdir -p %s' %os.path.dirname(params.temp))
        shell('samtools view -h -q {params.mapq} {input} | samtools view -bT {hpv_ref} -o {params.temp}')
        shell('samtools sort -o {output} -@ {threads} -T {params.pre} {params.temp}')
        shell('samtools index {output}; rm {params.temp}')

#--------------------------------------------------------------------------
# Note that TMAP automatically left aligns gaps and indels, so the GATK step
# from the original pipeline was removed.
#--------------------------------------------------------------------------

def amp_input(panel):   # check if we should use single panels or the universal panel amplicon bed
    if panel == 'universal':
        return list(config['panel'])
    else:
        return expand(config['amplicon_bed'] %'{hpvtype}', hpvtype=hpv_types)


rule amplicon_bed:
    input: amp_input
    output: 'refs/%s.amplicon.bed' %config['deliver_proj']
    run:
        shell('cat {input} > {output}')
        shell('sed -i "/^track/d" {output}') # remove lines that start with 'track'

rule length_bed:
    input: expand(config['len_bed'] %'{hpvtype}', hpvtype=hpv_types)
    output: 'refs/%s.len.bed' %config['deliver_proj']
    run:
        shell('cat {input} > {output}')

rule type_fastas:
    input: config['hpv_ref_nobreak']
    output: expand('refs/HPV{hpvtype}.fasta', hpvtype=hpv_types)
    run:
        for hpv in hpv_types:
            shell('grep -A1 "HPV%s_" {input} > refs/HPV%s.fasta' %(hpv, hpv))

rule variant_call:
    input: 
        rules.mapq_filter.output,
        rules.amplicon_bed.output,
        rules.length_bed.output
    output:
        'tvc/{sampleID}/TSVC_variants.vcf',
        'tvc/{sampleID}/{sampleID}.ptrim.bam'
    threads: 2
    params:
        pipe = config["vc_pipe"],
        out = ' tvc/{sampleID}',
        param = config["vc_param"],
        vc_bin = config["vc_bin"],
    run:
        shell('python {params.pipe} \
        --input-bam {input[0]} \
        --postprocessed-bam {output[1]} \
        --primer-trim-bed {input[1]} \
        --reference-fasta {hpv_ref} \
        --num-threads {threads} \
        --output-dir {params.out} \
        --parameters-file {params.param} \
        --bin-dir {params.vc_bin} \
        --region-bed {input[2]}')


rule adjust_padding:
    input: rules.variant_call.output[0]
    output: 'tvc_vcf/{sampleID}.tvc_no_pad.vcf'
    params: temp = '{sampleID}.temp.vcf'
    run:
        if config['padding'] == False:
            shell('cd tvc_vcf; ln -s ../{input} {wildcards.sampleID}.tvc_no_pad.vcf')
        else:
            vcf = open(input[0], 'r')
            outfile = open(output[0], 'w')
            need_sort = False
    
            for line in vcf:
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    type = line.split()[0]
                    hpv_len = bed[type]
                    loc = line.split()[1]
                    if int(loc) > hpv_len:
                        new_loc = int(loc) - hpv_len
                        outfile.write(line.replace(loc, str(new_loc), 1))
                        need_sort = True
                    else:
                        outfile.write(line)
            vcf.close()
            outfile.close()
    
            if need_sort == True:
                shell('vcf-sort -c {output} > {params.temp}')
                shell('mv {params.temp} {output}')

#--------------------------------------------------------------------------
rule hpv_bam:   # removes human reads
    input: 
        'tvc/{sampleID}/{sampleID}.ptrim.bam',
        rules.length_bed.output
    output: 'ptrim_hpv/{sampleID}.hpv.bam'
    run:
        shell('samtools view -h -L {input[1]} {input[0]} | samtools view -bS -o {output}')
        shell('samtools index {output}')

#--------------------------------------------------------------------------
rule pileup: # eventually update pipeline to use mpileup -aa to output all reference positions
    input: 
        rules.hpv_bam.output,
        rules.length_bed.output
    output: 'pileup/{sampleID}.pileup'
    run:
        shell('samtools mpileup -f {hpv_ref} -l {input[1]} {input[0]} > {output}')

#--------------------------------------------------------------------------
# Someday update this rule to use pyfaidx FastaVariant class
rule fasta:
    input:
        'pileup/{sampleID}.pileup', # pileup output
        'tvc_vcf/{sampleID}.tvc_no_pad.vcf', # adjust padding output
        expand('refs/HPV{hpv_type}.fasta', hpv_type=hpv_types)
    output:
        expand('fasta/{{sampleID}}_HPV{hpv_type}.fasta', hpv_type=hpv_types) # multiple type files for one sampleID
    run:
        # note vcf header is 0 because it comes after we've skipped 70 lines (as opposed to being the actual line 71)
        dvcf = pandas.read_table(input[1], skiprows=70, header=0)    
        dpile = pandas.read_table(input[0], names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
        types = list(set(dvcf['#CHROM'].tolist() + dpile['chrom'].tolist()))

        # create a fasta file for each HPV type in the project (so that the rule completes)
        # empty sample/types will be removed during the cat
        for hpv in hpv_types:
            print(hpv)
            fa = 'refs/HPV%s.fasta' %hpv
            fa_handle = open(fa, 'r')

            #record = SeqIO.parse(fa_handle, 'fasta').next()   # just read first record
            # TODO - SeqIOparse().next is throwing an error about no "next" attribute.
            seq = ''
            for record in SeqIO.parse(fa_handle, 'fasta'):
                seq = str(record.seq)
                break # this also takes just the first record, it's just longer
            fa_handle.close()
            
            # now start looking for SNPs and deletions as per original pipeline
            dt = dvcf[dvcf['#CHROM'] == 'HPV%s_Ref' %hpv].copy()
            
            if config['padding'] == True:
                newseq = seq[:len(seq)-400] # start with the ref sequence and add SNPs
            else:
                newseq = str(seq)
            for idx, row in dt.iterrows():
                (pos, ref, alt) = (int(row['POS']), row['REF'], row['ALT'].split(',')[0])
                # Look for SNPs
                # TODO:  add test for TVC REF matching the REF in the seq string            
                if len(ref) == len(alt):
                    counter = 0
                    while counter < len(ref):
                        newseq = newseq[:pos-1+counter] + alt[counter] + newseq[pos+counter:]
                        counter += 1

                # Look for deletions
                elif (len(ref) > len(alt)) and (len(alt) > 1):
                    counter = 1 # the first nt is the same as the reference, so start at +1
                    while counter < len(ref):
                        newseq = newseq[:pos-1] + '-' + newseq[pos:]
                        counter += 1

                # Skip insertions and other types of variation
                else:
                    continue

            # now check pileup to make sure there was enough coverage at each location
            dp = dpile[dpile['chrom'] == 'HPV%s_Ref' %hpv].copy()

            # Account for zero coverage in pileup (position is missing) by adding zeros
            dp.set_index('pos', inplace=True)
            allpos = list(range(len(seq)+1))
            allpos.pop(0)
            dp = dp.reindex(allpos).fillna(0) 
            dp['cov'] = dp['cov'].astype(int)

            print(newseq)

            # take padding into account if necessary
            if config['padding'] == True:
                dp['adj_pos'] = dp.index
                dp['adj_pos'] = dp['adj_pos'].apply(lambda x: (int(x) - (len(seq)-400)) if x > (len(seq)-400) else int(x))
                x = dp.groupby('adj_pos')['cov'].sum()
                
                # base calls with less than min_depth are called as N
                x = x[x < config['min_read']]
                print(len(x))
                for pos, depth in x.iteritems():
                    newseq = newseq[:pos-1] + 'N' + newseq[pos:]
            
            else:
                x = dp[dp['cov'] < 4].copy()
                for pos, depth in x.iterrows():
                    newseq = newseq[:pos-1] + 'N' + newseq[pos:]
            
            print(newseq)

            # output a fasta file for each type in each sample
            outfile = open('fasta/%s_HPV%s.fasta' %(wildcards.sampleID, hpv), 'w')
            outfile.write('>%s_HPV%s\n' %(wildcards.sampleID, hpv))
            outfile.write(newseq + '\n')
            outfile.close()


#--------------------------------------------------------------------------
rule fasta_cat:
    input: expand('fasta/{sampleID}_HPV{{hpv_type}}.fasta', sampleID=sampleIDs)
    output: 'cat_fasta/%s.HPV{hpv_type}.fasta' %config['deliver_proj']
    run:
        shell('cat {input} > {output}')


rule fasta_n:
    input: rules.fasta_cat.output
    output: 'reports/fasta/%s.HPV{hpv_type}.N-%d.fasta' %(config['deliver_proj'], config['fasta_n'])
    run:
        maxN = int(config['fasta_n'])
        keep = []
        for record in SeqIO.parse(input[0], 'fasta'):
            d = collections.Counter(record.seq) # creates a dictionary w/counts for each character found
            if 'N' in d.keys():
                if (float(d['N'])/len(record.seq) * 100) < maxN:
                    keep.append(record)
            else:
                keep.append(record)

        # SeqIO.write doesn't let you set wrap width, so use FastaIO directly
        outfile = open(output[0], 'w')
        fasta_out = FastaIO.FastaWriter(outfile, wrap=None)
        fasta_out.write_file(keep)
        outfile.close()



#--------------------------------------------------------------------------
rule type_summary:
    input: expand('pileup/{sampleID}.pileup', sampleID=d.keys())
    output: 'reports/type_summary.tsv'
    run:
        stacks = []
        # iterate through each sample pileup and pull out which types were found
        for sample in input:
            df = pandas.read_table(sample, names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
            df['sampleID'] = sample.split('/')[-1].split('.')[0]
            # note this does not consider padding yet!
            df = df[df['cov'] >= int(config['min_read'])]
            x = df.groupby(['sampleID', 'chrom'])['pos'].count()
            y = x.unstack()
            stacks.append(y) 

        dfs = pandas.concat(stacks).fillna(0)
        dfs.to_csv(output[0], sep='\t')


