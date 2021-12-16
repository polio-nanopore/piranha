import os
import collections
from Bio import SeqIO
import yaml
import json

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
import piranha.analysis.get_haplotypes

rule all:
    input:
        # expand(os.path.join(config[KEY_OUTDIR],"{taxid}","haplotype","pseudo_aln.fasta"), taxid=config["taxa_present"])
        expand(os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","variants.vcf"), taxid=config["taxa_present"]),
        os.path.join(config[KEY_OUTDIR],"consensus_sequences","cns.prompt.txt"),
        os.path.join(config[KEY_OUTDIR],"variation_info.json")

rule files:
    params:
        ref=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fasta"),
        reads=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fastq")

rule minimap2:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref,
    output:
        sam = os.path.join(config[KEY_OUTDIR],"{taxid}","polishing","mapped.ref.sam")
    shell:
        """
        minimap2 -ax map-ont {input.ref:q} {input.reads:q} > {output.sam:q}
        """
rule sort_index:
    input:
        sam = rules.minimap2.output.sam
    output:
        bam = os.path.join(config[KEY_OUTDIR],"{taxid}","polishing","mapped.sorted.bam"),
        index = os.path.join(config[KEY_OUTDIR],"{taxid}","polishing","mapped.sorted.bam.bai")
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
    params:
        outdir=os.path.join(config[KEY_OUTDIR],"{taxid}","medaka"),
        cns_mod = os.path.join(config[KEY_OUTDIR],"{taxid}","polishing","racon.mod.fasta")
    output:
        probs = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.sam:q} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft:q} > {params.cns_mod:q}
            medaka_consensus -i {input.basecalls:q} -d {params.cns_mod:q} -o {params.outdir:q} -t 2
        else
            touch {output.consensus:q}
        fi
        """

rule medaka_variant:
    input:
        probs = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus_probs.hdf"),
        ref = rules.files.params.ref,
        bam = rules.sort_index.output.bam
    params:
        temp_vcf = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","variants.pre.vcf")
    output:
        vcf = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","variants.vcf")
    shell:
        """
        medaka variant {input.ref:q} {input.probs:q} {params.temp_vcf:q} &&
        medaka tools annotate {params.temp_vcf:q} {input.ref:q} {input.bam:q} {output[0]}
        """

rule sam_to_seq:
    input:
        sam = rules.minimap2.output.sam,
        ref = rules.files.params.ref
    output:
        fasta = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","pseudoaln.fasta")
    shell:
        """
        gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q}
        """

rule get_variant_data:
    input:
        expand(rules.files.params.reads, taxid=config["taxa_present"]),
        expand(rules.sam_to_seq.output.fasta,taxid=config["taxa_present"])
    output:
        json = os.path.join(config[KEY_OUTDIR],"variation_info.json")
    run:
        variation_dict = {}
        for taxon in config["taxa_present"]:
            ref = os.path.join(config[KEY_TAXA_OUTDIR],f"{taxon}.fasta")
            fasta = os.path.join(config[KEY_OUTDIR],f"{taxon}","medaka","pseudoaln.fasta")
            var_dict = get_haplotypes.get_variation_pcent(ref,fasta)
            variation_dict[taxon] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))

rule get_haplotype_data:
    input:
        ref = rules.files.params.ref,
        fasta = rules.sam_to_seq.output.fasta,
        vcf = rules.medaka_variant.output.vcf,
        reads = rules.files.params.reads
    output:
        haplotypes = os.path.join(config[KEY_OUTDIR],"{taxid}","haplotypes.csv")
    run:
        get_haplotypes.get_haplotypes(input.ref,input.fasta,input.vcf,input.reads,
                        output,outdir,config[KEY_MIN_READS])
        
        
rule gather:
    input:
        expand(os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","variants.vcf"), taxid=config["taxa_present"])
    output:
        os.path.join(config[KEY_OUTDIR],"consensus_sequences","cns.prompt.txt")
    shell:
        """
        touch {output[0]}"""