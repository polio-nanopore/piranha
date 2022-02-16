import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.analysis.parse_paf import parse_paf_file
from piranha.analysis.filter_lengths import filter_reads_by_length
##### Target rules #####

rule all:
    input:
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","haplotypes.csv"), barcode=config["barcodes"]),
        # expand(os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html"), barcode=config["barcodes"]),
        # expand(os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences","cns.prompt.txt"), barcode=config["barcodes"]),
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html"), barcode=config["barcodes"]),
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","consensus_sequences.fasta"), barcode=config["barcodes"])


rule assess_haplotypes:
    input:
        snakefile = os.path.join(workflow.current_basedir,"haplotype_reader.smk"),
        taxa = rules.assess_broad_diversity.output.taxa,
        yaml = rules.assess_broad_diversity.output.new_config
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    threads: workflow.cores*0.5
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}_haplotype.smk.log")
    output:
        var_data = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","variation_info.json"),
        csv = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","haplotypes.csv"),
        config = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","haplotype_config.yaml")
    run:
        print(green(f"Assessing haplotypes for {params.barcode}"))
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "{config[log_string]} "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} tempdir={params.tempdir:q} "
                    "--cores {threads} &> {log:q}")

rule haplotype_consensus:
    input:
        yaml = rules.assess_haplotypes.output.config,
        snakefile = os.path.join(workflow.current_basedir,"haplotype_consensus.smk")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    output:
        fasta = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_data","consensus_sequences.fasta")
    run:
        print(green("Running snipit pipeline."))

        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                "--forceall "
                "{config[log_string]} "
                f"--directory '{config[KEY_TEMPDIR]}' "
                "--configfile {input.yaml:q} "
                "--config outdir={params.outdir:q}  tempdir={params.tempdir:q} "
                "--cores {workflow.cores} ")


rule generate_snipit:
    input:
        fasta = rules.haplotype_consensus.output.fasta,
        yaml = rules.assess_haplotypes.output.config,
        snakefile = os.path.join(workflow.current_basedir,"snipit.smk")
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}"),
        tempdir = os.path.join(config[KEY_TEMPDIR],"{barcode}")
    output:
        txt = os.path.join(config[KEY_OUTDIR],"{barcode}","snipit.txt")
    run:
        print(green("Running snipit pipeline."))
        print(input.fasta)
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                "--forceall "
                "{config[log_string]} "
                f"--directory '{config[KEY_TEMPDIR]}' "
                "--configfile {input.yaml:q} "
                "--config fasta={input.fasta:q} outdir={params.outdir:q}  tempdir={params.tempdir:q} "
                "--cores {workflow.cores} ")

rule generate_report:
    input:
        haplotype_info = rules.assess_haplotypes.output.csv,
        read_lengths = rules.filter_by_length.output.lengths,
        variation_info = rules.assess_haplotypes.output.var_data,
        yaml = rules.assess_haplotypes.output.config,
        snipit = rules.generate_snipit.output.txt
    params:
        barcode = "{barcode}",
    output:
        html = os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html")
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)

        make_report(output.html,input.read_lengths,input.variation_info,input.haplotype_info,params.barcode,config_loaded)


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
