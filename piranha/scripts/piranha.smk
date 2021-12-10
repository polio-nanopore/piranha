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
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample","refs_present.csv"), barcode=config["barcodes"]),
        expand(os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample","consensus_sequences","cns.prompt.txt"), barcode=config["barcodes"]),
        # expand(os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html"), barcode=config["barcodes"]),
        # expand(os.path.join(config[KEY_OUTDIR],"{barcode}","consensus_sequences","cns.prompt.txt"), barcode=config["barcodes"])


rule bin_files:
    input:
        fastq = os.path.join(config[KEY_READDIR], "{barcode}", f"fastq_runid_{config[KEY_RUNID]}_0.fastq")
    params:
        file_path = os.path.join(config[KEY_READDIR], "{barcode}")
    output:
        binned = os.path.join(config[KEY_OUTDIR],"{barcode}","unfiltered","nanopore_reads.fastq")
    shell:
        """
        cat {params.file_path}/fastq_*.fastq > {output[0]}
        """

rule filter_by_length:
    input:
        os.path.join(config[KEY_OUTDIR],"{barcode}","unfiltered","nanopore_reads.fastq")
    output:
        fastq = os.path.join(config[KEY_OUTDIR],"{barcode}","read_length_filtered","filtered_reads.fastq"),
        lengths = os.path.join(config[KEY_OUTDIR],"{barcode}","read_length_filtered","lengths.txt")
    run:
        filter_reads_by_length(input[0],output.fastq,output.lengths,config)

rule map_reads:
    input:
        ref = config[KEY_REFERENCE_SEQUENCES],
        fastq = rules.filter_by_length.output.fastq
    params:
        barcode = "{barcode}"
    threads: workflow.cores
    output:
        paf = os.path.join(config[KEY_OUTDIR],"{barcode}","read_length_filtered","filtered_reads.paf")
    shell:
        """
        minimap2 -t {threads} -x map-ont --secondary=no --paf-no-hit \
        {input.ref:q} \
        {input.fastq:q} > {output:q}
        """

rule assess_sample_composition:
    input:
        map_file = rules.map_reads.output.paf,
        fastq = rules.filter_by_length.output.fastq,
        ref = config[KEY_REFERENCE_SEQUENCES]
    params:
        tax_out = os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample"),
        barcode = "{barcode}"
    output:
        new_config = os.path.join(config[KEY_OUTDIR],"{barcode}","config.yaml"),
        taxa = os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample","refs_present.csv")
    run:
        parse_paf_file(input.map_file,
                        input.fastq,
                        output.taxa,
                        config[KEY_REFERENCE_SEQUENCES],
                        params.tax_out,
                        output.new_config,
                        params.barcode,
                        config)

rule get_consensus_sequences:
    input:
        snakefile = os.path.join(workflow.current_basedir,"cns_runner.smk"),
        taxa = rules.assess_sample_composition.output.taxa,
        yaml = rules.assess_sample_composition.output.new_config
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample","consensus_sequences")
    threads: 4
    output:
        taxa = os.path.join(config[KEY_OUTDIR],"{barcode}","categorised_sample","consensus_sequences","cns.prompt.txt")
    run:
        print("Running consensus pipeline.")
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                    "--forceall "
                    "--configfile {input.yaml:q} "
                    "--config barcode={params.barcode} outdir={params.outdir:q} "
                    "--cores 4 && touch {output.taxa}")

# rule qc_cns_output:
#     input:
#     output:
#     run:

# rule generate_report:
#     input:
#         all_reads = os.path.join(config[KEY_OUTDIR],"{barcode}","unfiltered","nanopore_reads.fastq"),
#         classified_reads = os.path.join(config[KEY_OUTDIR],"{barcode}","classified","classified_reads.unfiltered.fastq"),
#         filtered_reads = os.path.join(config[KEY_OUTDIR],"{barcode}","read_length_filtered","filtered_reads.fastq"),
#         classified_filtered = os.path.join(config[KEY_OUTDIR],"{barcode}","classified","classified_reads.filtered.fastq"),
#         taxa = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_taxa","taxa_present.csv"),
#         prompt = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt"),
#         yaml = os.path.join(config[KEY_OUTDIR],"{barcode}","config.yaml"),
#     params:
#         barcode = "{barcode}",
#         cns_path = os.path.join(config[KEY_OUTDIR],"{barcode}","processed_taxa","consensus_sequences")
#     output:
#         html = os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html")
#     run:
#         # for seq in os.path.join(config[KEY_OUTDIR],"{barcode}","processed_taxa","consensus_sequences"):
#         # if file empty, remove, else add to table
#         data_for_report = {}
#         report.get_read_counts(input.all_reads, input.classified_reads, input.filtered_reads, input.classified_filtered,data_for_report)
#         report.load_histograms_svg(input.histogram1, input.histogram2, data_for_report)

#         to_remove = report.data_for_table(input.taxa,params.cns_path, data_for_report)

#         with open(input.yaml, 'r') as f:
#             config_loaded = yaml.safe_load(f)

#         report.make_report(output.html,config_loaded,data_for_report,params.barcode)
#         # for i in to_remove:
#         #     shell(f"rm {i}")

# rule publish:
#     input:
#         report = os.path.join(config[KEY_OUTDIR],"{barcode}","analysis_report.html"),
#         fkrona = os.path.join(config[KEY_OUTDIR],"{barcode}","classified","krona.filtered.html"),
#         ukrona = os.path.join(config[KEY_OUTDIR],"{barcode}","classified","krona.unfiltered.html"),
#         cns= os.path.join(config[KEY_OUTDIR],"{barcode}","processed_taxa","consensus_sequences","cns.prompt.txt")
#     params:
#         barcode = "{barcode}",
#         pcns_path = os.path.join(config["publish_path"],"{barcode}","consensus_sequences")
#     output:
#         preport = os.path.join(config["publish_path"],"{barcode}","analysis_report.html"),
#         rreport =os.path.join(config["repo_path"],"{barcode}","analysis_report.html"),
#         pcns = os.path.join(config["publish_path"],"{barcode}","consensus_sequences","cns.prompt.txt"),
#         rcns = os.path.join(config["repo_path"],"{barcode}","consensus_sequences","cns.prompt.txt"),
#         pfkrona = os.path.join(config["publish_path"],"{barcode}","krona.filtered.html"),
#         pukrona = os.path.join(config["publish_path"],"{barcode}","krona.unfiltered.html"),
#         rfkrona = os.path.join(config["repo_path"],"{barcode}","krona.filtered.html"),
#         rukrona = os.path.join(config["repo_path"],"{barcode}","krona.unfiltered.html")
#     run:
#         shell(
#             """
#             cp {input.report} {output.preport}
#             cp {input.report} {output.rreport}
#             cp {input.fkrona} {output.pfkrona}
#             cp {input.fkrona} {output.rfkrona}
#             cp {input.ukrona} {output.pukrona}
#             cp {input.ukrona} {output.rukrona}
#             cp {input.cns} {output.pcns}
#             cp {input.cns} {output.rcns}
#             """)
#         try:
#             shell("""
#             cp {params.cns_path}/*.fasta {params.rcns_path}
#             cp {params.cns_path}/*.fasta {params.pcns_path}
#             """)
#         except:
#             print("No consensus files to copy for barcode {params.barcode}")
