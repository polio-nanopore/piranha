import os
import collections
from Bio import SeqIO
import yaml
import json

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.analysis.get_haplotypes import *

"""
output files
-------------
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.ref.sam"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.sorted.bam"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.sorted.bam.bai"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","consensus_probs.hdf"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","cns.mod.fasta"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","consensus.fasta"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","variants.pre.vcf"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","variants.vcf"
config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","pseudoaln.fasta"
config[KEY_TEMPDIR],"assess_haplotypes","haplotypes_{taxid}.csv"

config[KEY_OUTDIR],"processed_data","variation_info.json"
config[KEY_OUTDIR],"processed_data","haplotypes.csv"
config[KEY_OUTDIR],"processed_data","haplotype_config.yaml"
config[KEY_OUTDIR],"processed_data","{taxid}_{haplotype}.fastq"
"""

rule all:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","variants.vcf"), taxid=config["taxa_present"]),
        os.path.join(config[KEY_OUTDIR],"processed_data","variation_info.json"),
        os.path.join(config[KEY_OUTDIR],"processed_data","haplotype_config.yaml")

rule files:
    params:
        ref=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fasta"),
        reads=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fastq")

rule minimap2:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref,
    log: os.path.join(config[KEY_TEMPDIR],"logs","{taxid}.minimap2.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.ref.sam")
    shell:
        """
        minimap2 -ax map-ont {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """
rule sort_index:
    input:
        sam = rules.minimap2.output.sam
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.sorted.bam"),
        index = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","mapped.sorted.bam.bai")
    shell:
        """
        samtools view -bS -F 4 {input.sam:q} | samtools sort -o {output[0]:q} &&
        samtools index {output.bam:q} {output.index:q}
        """

rule medaka_consensus:
    input:
        basecalls=rules.files.params.reads,
        draft=rules.files.params.ref,
        sam= rules.minimap2.output.sam
    log: os.path.join(config[KEY_TEMPDIR],"logs","{taxid}.medaka_consensus.log")
    params:
        outdir=os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}"),
        cns_mod = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","cns.mod.fasta")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.sam:q} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft:q} > {params.cns_mod:q}
            medaka_consensus -i {input.basecalls:q} -d {params.cns_mod:q} -o {params.outdir:q} -t 2 &> {log:q}
        else
            touch {output.consensus:q}
        fi
        """

rule medaka_variant:
    input:
        probs = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","consensus_probs.hdf"),
        ref = rules.files.params.ref,
        bam = rules.sort_index.output.bam
    log: os.path.join(config[KEY_TEMPDIR],"logs","{taxid}.medaka_variant.log")
    params:
        temp_vcf = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","variants.pre.vcf")
    output:
        vcf = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","variants.vcf")
    shell:
        """
        medaka variant {input.ref:q} {input.probs:q} {params.temp_vcf:q}  &> {log} 
        medaka tools annotate {params.temp_vcf:q} {input.ref:q} {input.bam:q} {output[0]} &> {log}
        """

rule sam_to_seq:
    input:
        sam = rules.minimap2.output.sam,
        ref = rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{taxid}.gofasta.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","{taxid}","pseudoaln.fasta")
    shell:
        """
        gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q} &> {log}
        """

rule get_variant_data:
    input:
        expand(rules.files.params.reads, taxid=config["taxa_present"]),
        expand(rules.sam_to_seq.output.fasta,taxid=config["taxa_present"])
    output:
        json = os.path.join(config[KEY_OUTDIR],"processed_data","variation_info.json")
    run:
        # this is for making a figure
        variation_dict = {}
        for taxon in config["taxa_present"]:
            ref = os.path.join(config[KEY_TAXA_OUTDIR],f"{taxon}.fasta")
            fasta = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes",f"{taxon}","pseudoaln.fasta")
            
            var_dict = get_variation_pcent(ref,fasta)
            variation_dict[taxon] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))

rule get_haplotype_data:
    input:
        fasta = rules.sam_to_seq.output.fasta,
        vcf = rules.medaka_variant.output.vcf,
        reads = rules.files.params.reads
    params:
        taxon = "{taxid}",
        outdir = os.path.join(config[KEY_OUTDIR],"processed_data")
    output:
        haplotypes = os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","haplotypes_{taxid}.csv")
    run:
        get_haplotypes(input.fasta,input.vcf,input.reads,
                        output.haplotypes,params.outdir,params.taxon,config[KEY_MIN_READS],config[KEY_MIN_PCENT])

rule gather_haplotypes:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"assess_haplotypes","haplotypes_{taxid}.csv"), taxid=config["taxa_present"])
    output:
        csv = os.path.join(config[KEY_OUTDIR],"processed_data","haplotypes.csv"),
        config = os.path.join(config[KEY_OUTDIR],"processed_data","haplotype_config.yaml")
    run:
        gather_haplotype_data(input,output.csv,output.config,config)
