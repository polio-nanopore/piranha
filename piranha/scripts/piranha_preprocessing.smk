#!/usr/bin/env python3

import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *
from piranha.analysis.preprocessing import *
from piranha.analysis.filter_lengths import filter_reads_by_length
##### Target rules #####

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt"), barcode=config[KEY_BARCODES]),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config[KEY_BARCODES])

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
        paf = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","filtered_reads.paf")
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
        hits = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","hits_reads.csv")
    run:
        parse_paf_file(input.map_file,
                        output.csv,
                        output.hits,
                        config[KEY_REFERENCE_SEQUENCES],
                        params.barcode,
                        config)

rule write_hit_fastq:
    input:
        csv = rules.assess_broad_diversity.output.csv,
        hits = rules.assess_broad_diversity.output.hits,
        fastq = rules.filter_by_length.output.fastq
    params:
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt")
    run:
        shell("touch {output.txt:q}")
        print(green(params.barcode))
        to_write = write_out_fastqs(input.csv,input.hits,input.fastq,params.outdir,config)

        write_out_ref_fasta(to_write,config[KEY_REFERENCE_SEQUENCES],params.outdir)


rule gather_diversity_report:
    input:
        refs = expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config[KEY_BARCODES]),
        txt = expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt"), barcode=config[KEY_BARCODES])
    output:
        refs= os.path.join(config[KEY_TEMPDIR],SAMPLE_COMPOSITION),
        summary = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG)
    run:
        barcode_config = diversity_report(input.refs,output.refs,output.summary,config)
        
        with open(output.yaml, 'w') as fw:
            yaml.dump(barcode_config, fw) 
        print(yellow("-----------------------"))
