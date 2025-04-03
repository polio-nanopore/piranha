#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.stool_functions import *
from piranha.report.make_report import make_sample_report
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
##### Target rules #####
"""
input files
os.path.join(config[KEY_OUTDIR],PREPROCESSING_SUMMARY)
os.path.join(config[KEY_OUTDIR],"hits.csv")
os.path.join(config[KEY_OUTDIR],SAMPLE_COMPOSITION)
"""
rule all:
    input:
        os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS),
        expand(os.path.join(config[KEY_OUTDIR],"barcode_reports","{barcode}_report.html"), barcode=config[KEY_BARCODES]),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.fasta"), barcode=config[KEY_BARCODES])

rule files:
    params:
        composition=os.path.join(config[KEY_TEMPDIR],SAMPLE_COMPOSITION),
        summary=os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG)


rule estimate_haplotypes:
    input:
        snakefile = os.path.join(workflow.current_basedir,"piranha_haplotype.smk"),
        yaml = rules.files.params.yaml,
        prompt = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}"),
        publish_dir = os.path.join(config[KEY_OUTDIR],"published_data","{barcode}")
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_haplotype.smk.log")
    output:
        yaml = os.path.join(config[KEY_TEMPDIR],"{barcode}",HAPLOTYPING_CONFIG)
    run:
        if config[KEY_RUN_HAPLOTYPING]:
            sample = get_sample(config[KEY_BARCODES_CSV],params.barcode)
            print(green(f"Calculating haplotypes for {sample} ({params.barcode})"))
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "--forceall "
                        "--rerun-incomplete "
                        "{config[log_string]} "
                        "--configfile {input.yaml:q} "
                        "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                        f"sample='{sample}' "
                        "--cores {threads} &> {log:q}")
        else:
            shell(
                """
                cp {input.yaml:q} {output.yaml:q}
                """)

rule generate_consensus_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"piranha_consensus.smk"),
        yaml = rules.estimate_haplotypes.output.yaml
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_consensus.smk.log")
    output:
        yaml = os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_config.yaml")
    run:
        sample = get_sample(config[KEY_BARCODES_CSV],params.barcode)
        print(green(f"Calculating consensus sequences for {sample} ({params.barcode})"))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--rerun-incomplete "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    f"sample='{sample}' "
                    "--cores {threads} &> {log:q}")


rule curate_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"piranha_curate.smk"),
        yaml = rules.generate_consensus_sequences.output.yaml
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}"),
        variant_dir = os.path.join(config[KEY_TEMPDIR],"{barcode}","variant_calls"),
        publish_dir = os.path.join(config[KEY_OUTDIR],"published_data","{barcode}")
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_curate.smk.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.fasta"),
        csv= os.path.join(config[KEY_TEMPDIR],"{barcode}","variants.csv"),
        masked =  os.path.join(config[KEY_TEMPDIR],"{barcode}","masked_variants.csv")
    run:
        sample = get_sample(config[KEY_BARCODES_CSV],params.barcode)
        print(green(f"Curating consensus sequences & variants for {sample} ({params.barcode})"))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--rerun-incomplete "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    f"sample='{sample}' "
                    "--cores {threads} &> {log:q}")
        shell(
            """
            cp {params.variant_dir}/*.vcf {params.publish_dir}
            """)

rule generate_variation_info:
    input:
        snakefile = os.path.join(workflow.current_basedir,"piranha_variation.smk"),
        yaml = rules.generate_consensus_sequences.output.yaml,
        variants = rules.curate_sequences.output.csv,
        fasta = os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.fasta")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_variation.smk.log")
    output:
        json = os.path.join(config[KEY_TEMPDIR],"{barcode}","variation_info.json"),
        json_mask = os.path.join(config[KEY_TEMPDIR],"{barcode}","mask_info.json")
    run:
        # decide if we want 1 per haplotyde or 1 per ref group, will need mods either way
        sample = get_sample(config[KEY_BARCODES_CSV],params.barcode)
        print(green(f"Gathering variation info for {sample} ({params.barcode})"))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--rerun-incomplete "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    f"sample='{sample}' "
                    "--cores {threads} &> {log:q}")


rule mask_consensus_sequences:
    input:
        mask_json = rules.generate_variation_info.output.json_mask,
        fasta = rules.curate_sequences.output.fasta
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.masked.fasta")
    run:
        mask_low_coverage(input.mask_json, input.fasta,output.fasta)

rule gather_consensus_sequences:
    input:
        composition = rules.files.params.composition,
        fasta = expand(rules.mask_consensus_sequences.output.fasta, barcode=config[KEY_BARCODES])
    params:
        publish_dir = os.path.join(config[KEY_OUTDIR],"published_data")
    output:
        fasta = os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS)
    run:
        print(green("Gathering fasta files"))
        # needs mod for checking & merging identical cns within ref group
        # also header now needs hap parsing & hap->CNS mapping
        gather_fasta_files(input.composition, config[KEY_BARCODES_CSV], input.fasta,config[KEY_ALL_METADATA], output[0],params.publish_dir,config)

rule generate_report:
    input:
        consensus_seqs = rules.gather_consensus_sequences.output.fasta,
        variation_info = rules.generate_variation_info.output.json,
        masked_variants = rules.curate_sequences.output.masked,
        variants = rules.curate_sequences.output.csv,
        cns_yaml = rules.generate_consensus_sequences.output.yaml,
        yaml = rules.files.params.yaml
    params:
        outdir = os.path.join(config[KEY_OUTDIR],"barcode_reports"),
        barcode = "{barcode}",
    output:
        html = os.path.join(config[KEY_OUTDIR],"barcode_reports","{barcode}_report.html")
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)
        with open(input.cns_yaml, 'r') as f:
            cns_config_loaded = yaml.safe_load(f)

        #var dict now has total on it- so can infer from var dict which sites have been masked out
        make_sample_report(output.html,
                            input.variation_info,
                            input.consensus_seqs,
                            input.masked_variants,
                            params.barcode,
                            cns_config_loaded,
                            config_loaded)


