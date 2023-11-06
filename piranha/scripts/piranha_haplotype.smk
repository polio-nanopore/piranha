

import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.haplotyping_functions import *
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

rule minimap_for_bam:
    input:
    output:
    shell:

rule freebayes:
    input:
    output:
    shell:

rule flopp:
    input:
    output:
    shell:

rule merge_close_haplo:
    input:
    output:
    shell:

rule haplo_qc:
    input:
    output:
    shell:

rule generate_partition_files:
    input:
    output:
    shell:


rule medaka_haploid_variant:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    params:
        model = config[KEY_MEDAKA_MODEL],
        outdir =  os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus_probs.hdf"),
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","medaka.vcf"),
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus.fasta"),
        pub_vcf = os.path.join(config[KEY_TEMPDIR],"variant_calls","{reference}.vcf")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.hapoid_variant.log")
    shell:
        """
        [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}
        if [ -s {input.ref:q} ]
        then
            medaka_haploid_variant -i {input.reads:q} \
                                -r {input.ref:q} \
                                -o {params.outdir:q} \
                                -f -x && \
            medaka stitch {output.probs:q} {input.ref:q} {output.cns:q}
        else
            touch {output.cns:q}
            touch {output.probs:q}
            touch {output.vcf:q}
        fi
        cp {output.vcf:q} {output.pub_vcf:q}
        """

rule medaka_haploid_variant_cns:
    input:
        reads=rules.files.params.reads,
        ref=rules.medaka_haploid_variant.output.cns
    params:
        model = config[KEY_MEDAKA_MODEL],
        outdir = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus_probs.hdf"),
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","medaka.vcf"),
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus.fasta")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}_cns.hapoid_variant.log")
    shell:
        """
        [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}
        if [ -s {input.ref:q} ]
        then
            medaka_haploid_variant -i {input.reads:q} \
                                -r {input.ref:q} \
                                -o {params.outdir} \
                                -f -x && \
            medaka stitch {output.probs} {input.ref} {output.cns}
        else
            touch {output.cns:q}
            touch {output.probs:q}
            touch {output.vcf:q}
        fi
        """
