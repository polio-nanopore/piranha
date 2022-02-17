#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.stool_functions import *
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
        os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta"),
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences.fasta"), barcode=config[KEY_BARCODES])

rule files:
    params:
        composition=os.path.join(config[KEY_OUTDIR],SAMPLE_COMPOSITION),
        summary=os.path.join(config[KEY_OUTDIR],PREPROCESSING_SUMMARY)


rule generate_consensus_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"consensus.smk"),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG),
        prompt = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    threads: workflow.cores*0.5
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_consensus.smk.log")
    output:
        os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences.fasta"),
        os.path.join(config[KEY_OUTDIR],"{barcode}","variants.csv")
    run:
        print(green(f"Generating consensus sequences for {params.barcode}"))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    "--cores {threads} &> {log:q}")

rule gather_consensus_sequences:
    input:
        composition = rules.files.params.composition,
        fasta = expand(os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences.fasta"), barcode=config[KEY_BARCODES])
    output:
        fasta = os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta")
    run:
        gather_fasta_files(input.composition, config[KEY_BARCODES_CSV], input.fasta, output[0])

rule generate_snipits:
    input:
        snakefile = os.path.join(workflow.current_basedir,"snipit.smk"),
        cns = os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences.fasta"),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG)
    params:
        barcode = "{barcode}"
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"{barcode}","snipit.txt")
    run:
        print(green("Running snipit pipeline."))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                "--forceall "
                "{config[log_string]} "
                f"--directory '{config[KEY_TEMPDIR]}' "
                "--configfile {input.yaml:q} "
                "--config cns={input.cns:q} outdir={params.outdir:q}  tempdir={params.tempdir:q} "
                "--cores {workflow.cores} ")

# rule generate_report:
#     input:
#         haplotype_info = rules.assess_haplotypes.output.csv,
#         read_lengths = rules.filter_by_length.output.lengths,
#         variation_info = rules.assess_haplotypes.output.var_data,
#         yaml = rules.assess_haplotypes.output.config,
#         snipit = rules.generate_snipit.output.txt
#     params:
#         barcode = "{barcode}",
#     output:
#         html = os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html")
#     run:
#         with open(input.yaml, 'r') as f:
#             config_loaded = yaml.safe_load(f)

#         make_report(output.html,input.read_lengths,input.variation_info,input.haplotype_info,params.barcode,config_loaded)


# rule get_overall_report:
#     input:
#         yaml = rules.assess_haplotypes.output.config
#     params:
#         barcode = "{barcode}",
#         outdir = os.path.join(config[KEY_OUTDIR],"{barcode}")
#     threads: workflow.cores*0.5
#     log: os.path.join(config[KEY_TEMPDIR],"{barcode}_consensus.smk.log")
#     output:
#         html = os.path.join(config[KEY_OUTDIR],"piranha_report.html")
#     run:
#         with open(input.yaml, 'r') as f:
#             config_loaded = yaml.safe_load(f)
#         print(green(f"Generating haplotype consensus sequences for {params.barcode}."))
#         for h in config_loaded["haplotypes"]:
#             print(f"- {h}")
#         print("----------------")
#         shell("snakemake --nolock --snakefile {input.snakefile:q} "
#                     "--forceall "
#                     "{config[log_string]} "
#                     "--configfile {input.yaml:q} "
#                     "--config barcode={params.barcode} outdir={params.outdir:q} "
#                     "--cores {threads} && touch {output.taxa}")

