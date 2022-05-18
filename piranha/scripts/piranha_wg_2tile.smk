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

rule all:
    input:
        os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS),
        expand(os.path.join(config[KEY_OUTDIR],"barcode_reports","{barcode}_report.html"), barcode=config[KEY_BARCODES]),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.fasta"), barcode=config[KEY_BARCODES])


rule files:
    params:
        composition=os.path.join(config[KEY_TEMPDIR],SAMPLE_COMPOSITION),
        summary=os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY)


rule generate_consensus_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"consensus.smk"),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG),
        prompt = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_consensus.smk.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"{barcode}","consensus_sequences.fasta"),
        csv= os.path.join(config[KEY_TEMPDIR],"{barcode}","variants.csv"),
        masked =  os.path.join(config[KEY_TEMPDIR],"{barcode}","masked_variants.csv"),
        json = os.path.join(config[KEY_TEMPDIR],"{barcode}","variation_info.json")
    run:
        print(green(f"Generating consensus sequences for {params.barcode}"))
        sample = get_sample(config[KEY_BARCODES_CSV],params.barcode)
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--rerun-incomplete "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    f"sample='{sample}' "
                    "--cores {threads} &> {log:q}")