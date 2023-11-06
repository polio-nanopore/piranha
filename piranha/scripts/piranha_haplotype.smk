

import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
import piranha.analysis.haplotyping_functions as hf
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta"),
        os.path.join(config[KEY_TEMPDIR],"variants.csv"),
        os.path.join(config[KEY_TEMPDIR],"masked_variants.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"variant_calls","{reference}.vcf"), reference=REFERENCES),
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg"), reference=REFERENCES),

rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule rasusa:
    input:
        reads= rules.files.params.reads,
        ref=rules.files.params.ref
    params:
        depth = config[KEY_HAPLOTYPE_SAMPLE_SIZE]
    output:
        fastq= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","downsample.fastq")
    script:
        ref_len = 0
        for record in SeqIO.parse(input.ref,"fasta"):
            ref_len = len(record)
        shell("rasusa -i {input.reads:q} -c {params.depth:q} " + f"-g {ref_len}b" + " -o {output.reads_ds:q}")


rule minimap_for_bam:
    input:
        ref=rules.files.params.ref,
        fastq = rules.rasusa.output.fastq
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2_for_bam.log")
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","downsample.bam")
    shell:
        """
        minimap2 -t {threads} -ax asm20 --secondary=no  \
        {input.ref:q} \
        {input.fastq:q} | samtools sort -@ {threads} -O bam -o {output:q} &> {log:q}
        """

rule freebayes:
    input:
        bam = rules.minimap_for_bam.output.bam,
        ref = rules.files.params.ref
    output:
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","freebayes.vcf"),
        int_vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","freebayes_int.vcf")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.freebayes.log")
    params:
        allele_freq = config[KEY_MIN_ALLELE_FREQUENCY]
    run:
        shell("""
        freebayes -f {input.ref:q} \
                -F {params.alleleFreq} \
                --pooled-continuous \
                --no-indels \
                --no-mnps \
                --no-complex \
                {input.bam:q} > {output.int_vcf:q} 2> {log:q}
        """)
        
        hf.write_contig_headers_vcf(output.int_vcf,output.vcf)

rule flopp:
    input:
        bam = rules.minimap_for_bam.output.bam,
        vcf = rules.freebayes.output.vcf
    threads: workflow.cores
    params:
        ploidy = config[KEY_MAX_HAPLOTYPES],
        os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_partitions")
    output:
        flopp = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_output.flopp"),
        partition = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_partitions","partition.txt")
    shell:
        """
        VARIANTS=$(grep -v '#' {input.vcf} | wc -l)
        if [[ $VARIANTS -gt 0 ]]
        then
            flopp -b {input.bam:q} \
            -c {input.vcf:q} \
            -p {params.ploidy} -m \
            -o {output.flopp:q} \
            -t {threads} \
            -P {params.partitionPath} 
        else
            echo "No variants called - single reference haplotype"
            touch {output.flopp:q}
            touch {output.partition:q}
        fi
        """

rule merge_close_haplo:
    input:
        flopp = rules.flopp.output.flopp
    params:
        min_distance = config[KEY_MIN_HAPLOTYPE_DISTANCE]
    output:
        flopp = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","merged.flopp")
    run:
        ## merge haplotypes accross flopp file if within x snps

rule haplotype_qc:
    input:
        flopp = rules.flopp.output.flopp,
        partition = rules.flopp.output.partition,
        reads = rules.files.params.reads
    params:
        min_distance = config[KEY_MIN_HAPLOTYPE_DISTANCE],
        min_reads = config[KEY_MIN_HAPLOTYPE_DEPTH]
    output:
        haplotypes = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","curated_haplotypes.txt")
    run:
        ## qc steps
        """
        evenness statistic, standard deviation
        how many reads includued
        minimum reads for the haplotype to be used
        """

        #generate partition read files
        #need a file with list of haplotypes generated 
        #(i.e. name of read files produced)
