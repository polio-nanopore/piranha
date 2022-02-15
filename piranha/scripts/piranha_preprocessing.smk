import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.analysis.parse_paf import parse_paf_file,diversity_report
from piranha.analysis.filter_lengths import filter_reads_by_length
from piranha.report.make_report import make_report
##### Target rules #####

rule all:
    input:
        os.path.join(config[KEY_OUTDIR],"preprocessing_summary.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config["barcodes"])

rule gather_files:
    input:
        fastq = os.path.join(config[KEY_READDIR], "{barcode}", f"fastq_runid_{config[KEY_RUNID]}_0.fastq")
    params:
        file_path = os.path.join(config[KEY_READDIR], "{barcode}")
    output:
        fastq = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","nanopore_reads.fastq")
    shell:
        """
        cat {params.file_path}/fastq_*.fastq > {output[0]}
        """

rule filter_by_length:
    input:
        rules.gather_files.output.fastq
    output:
        fastq = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","filtered_reads.fastq")
    run:
        filter_reads_by_length(input[0],output.fastq,config)

rule map_reads:
    input:
        ref = config[KEY_REFERENCE_SEQUENCES],
        fastq = rules.filter_by_length.output.fastq
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}.minimap2_initial.log")
    output:
        paf = os.path.join(config[KEY_OUTDIR],"{barcode}","initial_processing","filtered_reads.paf")
    shell:
        """
        minimap2 -t {threads} -x map-ont --secondary=no --paf-no-hit \
        {input.ref:q} \
        {input.fastq:q} -o {output:q} &> {log:q}
        """

rule assess_broad_diversity:
    input:
        map_file = rules.map_reads.output.paf,
        ref = config[KEY_REFERENCE_SEQUENCES]
    params:
        barcode = "{barcode}"
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"),
        hits = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","hits_min_reads.csv")
    run:
        parse_paf_file(input.map_file,
                        output.csv,
                        output.hits,
                        config[KEY_REFERENCE_SEQUENCES],
                        params.barcode,
                        config)



rule gather_diversity_report:
    input:
        refs = expand(os.path.join(config[KEY_OUTDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config["barcodes"]),
        hits = expand(os.path.join(config[KEY_OUTDIR],"{barcode}","initial_processing","hits_min_reads.csv"), barcode=config["barcodes"])
    output:
        refs= os.path.join(config[KEY_OUTDIR],"sample_composition.csv"),
        summary = os.path.join(config[KEY_OUTDIR],"preprocessing_summary.csv"),
        hits = os.path.join(config[KEY_OUTDIR],"hits.csv")
    run:
        shell("cat {input.hits} > {output.hits:q}")

        refs_present = diversity_report(input.refs,output.refs,output.summary,config[KEY_BARCODES_CSV])

        print("----------------")
        print(green(f"Broad diversity assessed."))
        print("----------------")
