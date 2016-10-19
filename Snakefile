# vim: ft=python
import os
import glob
import itertools
import pysam

# --------------------------------------
# Original pipeline by Xijun Zhang
# Adapted for Snakemake by Mia Steinberg
# --------------------------------------

configfile: "hpv16_config.yaml"
workdir: os.environ['PWD']
shell.executable('bash')

#ANALYSIS_DIR = config["ANALYSIS_DIR"]

gatk = config["gatk"]
hg_fasta = config["hg_fasta"]
len_bed = config["len_bed"]
amp_bed = config["amp_bed"]
cov_dev = config["cov_dev"]

## ---- The parser will have to be customized for each run ---- ##
def parse_sampleID(filename):
    return filename.split('/')[-1].split('_')[2]

bamfiles = sorted(glob.glob(config['glob_path']), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(bamfiles, parse_sampleID):
    d[key] = list(value)

# We now have a dictionary with the sample ID as the keys and a list of 
# paths to the bam(s) as the value.
# e.g. d = {sample01: ['a/sample01.bam', 'b/sample01.bam'],
#           sample02: ['a/sample02.bam', 'c/sample02.bam'],
#           sample03: ['d/sample03.bam']}
         
# These rules run on the host node and are not submitted to the cluster.
localrules: all


#--------------------------------------------------------------------------
rule all:
    input:  expand('fasta/{sampleID}.fasta', sampleID=d.keys()), 
            expand('tvc_ann/{sampleID}.ann.vcf', sampleID=d.keys())


def link_rule_input_files(wildcards):
    return d[wildcards.sampleID]

#--------------------------------------------------------------------------
rule link:
    input: link_rule_input_files
    output: 'bams/{sampleID}.bam'
    run:
        if (len(input) > 1):
            tempbam = wildcards.sampleID + '.temp.bam'
            tempsam = wildcards.sampleID + '.temp.sam'

            shell('samtools merge {tempbam} {input}')
            
            # all SM fields in @RG must be identical
            # create a samfile with the fixed header
            sam = pysam.Samfile(tempbam, 'rb')
            header = sam.header
            for i in header['RG']:
                i['SM'] = wildcards.sampleID
            outfile = pysam.Samfile(tempsam, 'wh', header=header)
            
            # add the reads from the original merged bam
            shell('samtools view {tempbam} >> {tempsam}')
            
            # convert back to bam
            shell('samtools view -h {tempsam} > {output}')
            shell('rm {tempbam}; rm {tempsam}')
        else:
            shell('cd bams; ln -s {input} {wildcards.sampleID}.bam && touch -h {wildcards.sampleID}.bam')

#--------------------------------------------------------------------------
rule aq_filter:
    input: rules.link.output
    output: 'aq30/{sampleID}.filtered.bam'
    threads: 2
    params: temp = '{sampleID}.aq.temp'
    run:
        shell('samtools view -h -q 30 {input} | samtools view -bT {hg_fasta} -o {params.temp}')
        shell('samtools sort -o {output} -@ {threads} -T {wildcards.sampleID} {params.temp}')
        shell('samtools index {output}; rm {params.temp}')

#--------------------------------------------------------------------------
rule left_align:
    input: rules.aq_filter.output
    output: 'left/{sampleID}.left.bam'
    threads: 4
    run:
        shell('java -version')
        shell('java -Xmx8g -jar {gatk} -R {hg_fasta} -T LeftAlignIndels -I {input} -o {output}' )
        shell('samtools index {output}')
        shell('mv left/{wildcards.sampleID}.left.bai left/{wildcards.sampleID}.left.bam.bai')

#--------------------------------------------------------------------------
rule variant_call:
    input: rules.left_align.output
    output: 
        'tvc/{sampleID}/TSVC_variants.vcf',
        'tvc/{sampleID}/{sampleID}.ptrim.bam'
    params:        
        pipe = config["vc_pipe"],
        ptrim = 'tvc/{sampleID}/{sampleID}.ptrim.bam',
        out = ' tvc/{sampleID}',
        param = config["vc_param"],
        vc_bin = config["vc_bin"]
    threads: 4
    run:
        shell('python {params.pipe} \
        --input-bam {input} \
        --postprocessed-bam {params.ptrim} \
        --primer-trim-bed {amp_bed} \
        --reference-fasta {hg_fasta} \
        --output-dir {params.out} \
        --parameters-file {params.param} \
        --num-threads={threads} \
        --bin-dir {params.vc_bin} \
        --region-bed {len_bed}')

#--------------------------------------------------------------------------
rule hpv_bam:
    input: 'tvc/{sampleID}/{sampleID}.ptrim.bam'
    output: 'ptrim_hpv/{sampleID}.hpv.bam'
    run:
        shell('samtools view -h -L {len_bed} {input} | samtools view -bS -o {output}')
        shell('samtools index {output}')
       
#--------------------------------------------------------------------------
# GATK Haplotype caller is only used for compiling the fasta.
# It is not considered when making variant calls and will eventually be removed
# from the pipeline.
rule ht_caller:
    input: rules.hpv_bam.output
    output: 
        'hc_vcf/{sampleID}.hc.vcf',
        'hc_vcf_var/{sampleID}.hc.var.vcf'
    run:
        shell('java -Xmx10G -jar {gatk} \
                -T HaplotypeCaller \
                -R {hg_fasta} \
                -ploidy 1 \
                -I {input} \
                --emitRefConfidence BP_RESOLUTION \
                --output_mode EMIT_ALL_SITES \
                --allSitePLs \
                --dontUseSoftClippedBases \
                --variant_index_type LINEAR \
                --variant_index_parameter 128000 \
                -dt NONE \
                -L {len_bed} \
                -o {output[0]}')
        shell('java -Xmx8G -jar {gatk} \
                -T HaplotypeCaller \
                -R {hg_fasta} \
                -ploidy 1 \
                -I {input} \
                -stand_call_conf 30 \
                -stand_emit_conf 10 \
                --dontUseSoftClippedBases \
                -dt NONE \
                -L {len_bed} \
                -o {output[1]}')

#--------------------------------------------------------------------------
rule pileup:
    input: rules.hpv_bam.output
    output: 'pileup/{sampleID}.pileup'
    run:
        shell('samtools mpileup -f {hg_fasta} -l {len_bed} {input} > {output}')

#--------------------------------------------------------------------------
rule fasta:
    input: 
        tvc_file=rules.variant_call.output[0],
        hc_file=rules.ht_caller.output[0],
        pile_file=rules.pileup.output[0]
    output: 
        'fasta/{sampleID}.fasta',
        'gatk/{sampleID}_gatk.csv'
    run:
        outfile = open(output[0], 'w')
        flagfile = open(output[1], 'w')
        sample_seq = ""
        d = {}
        # create a dictionary with a key for each sequence location
        # the values are tuples of reference nt, alternate nt, and type of variation
        tvc = open(input.tvc_file, 'r')
        for line in tvc:
            if line.startswith("#"):
                continue
            vals = line.split('\t')
            (tvc_pos, tvc_ref, tvc_alt) = (int(vals[1]), vals[3], vals[4].split(',')[0])
            print(tvc_pos)
            # Look for SNPs
            if len(tvc_ref) == len(tvc_alt):
                counter = 0
                while counter < len(tvc_ref):
                    d[tvc_pos+counter] = (tvc_ref[counter], tvc_alt[counter], 'snp')
                    counter += 1

            # Look for deletions
            elif (len(tvc_ref) > len(tvc_alt)) and (len(tvc_alt) > 1):
                d[tvc_pos] = (tvc_ref[0], tvc_ref[0], 'del0')
                # first base in a deletion record, ref=alt
                counter = 1
                while counter < len(tvc_ref):
                    d[tvc_pos+counter] = (tvc_ref[counter], '-', 'del')
                    counter += 1

            else:
                # skip insertions and other types of variation
                print('skipping')
                continue
        tvc.close()
        
        hc = open(input.hc_file, 'r')
        for line in hc:
            if line.startswith('#'):
                continue
            vals = line.split('\t')
            zygosity = vals[9].split(':')[0].strip()
            (hc_pos, hc_ref, hc_alt) = (int(vals[1]), vals[3], vals[4].split(',')[0])
           
            pile = open(input.pile_file)
            dp = 0
            for line in pile:
                if int(line.split('\t')[1]) == hc_pos:
                    dp = int(line.split('\t')[3])
                    continue
            # take the base from the tvc file if variant is called
            if hc_pos in d.keys():
                sample_seq += d[hc_pos][1]
            # take the reference base from the hc file if read depth is > 6    
            elif (dp >= 6) and (zygosity == '0'):
                sample_seq += hc_ref[0]
            # make a note of the GATK SNPs that TVC did not call (and use ref)
            elif (dp >= 6) and (zygosity == '1'):
                flagfile.write('%s,%s,%s\n' %(wildcards.sampleID, hc_pos, hc_alt) )
                sample_seq  += hc_ref[0]
            # Use a N if the read depth is < 6    
            else:
                sample_seq += 'N'
        hc.close()

        outfile.write('>%s\n' %wildcards.sampleID)
        outfile.write(sample_seq + '\n')
        outfile.close()
        flagfile.close()

#--------------------------------------------------------------------------
rule annotate:
    input: rules.variant_call.output[0]
    output: 'tvc_ann/{sampleID}.ann.vcf'
    params: 
        snpeff = config["snpeff"],
        bed = config["snpeff_bed"],
        type = config["hpv_type"]
    run:
        shell('java -Xmx2g -jar {params.snpeff}/snpEff.jar \
                -ud 0 -interval {params.bed} \
                -c {params.snpeff}/snpEff.config {params.type} {input} > {output}')

#--------------------------------------------------------------------------
#onsuccess:
#    # run the qc scripts and copy deliverables to HPV Consortium
#    cov_dev = config["cov_dev"]
#    deliver = config["deliver_proj"]
#    shell('bash qc_report.sh {cov_dev} . {amp_bed} {hg_fasta} {deliver}')
