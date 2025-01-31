#!/usr/bin/env python3

import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *
from piranha.analysis.preprocessing import *
##### Target rules #####

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt"), barcode=config[KEY_BARCODES]),
        expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config[KEY_BARCODES])

rule filter_by_length:
    input:
    params:
        file_path = os.path.join(config[KEY_READDIR], "{barcode}"),
        barcode = "{barcode}"
    output:
        fastq = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","filtered_reads.fastq"),
        report = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","filter_report.csv")
    run:
        gather_filter_reads_by_length(params.file_path,params.barcode,output.fastq,output.report,config)

rule map_reads:
    input:
        ref = config[KEY_REFERENCE_SEQUENCES],
        fastq = rules.filter_by_length.output.fastq
    threads: workflow.cores
    params:
        minimap2_options = config[KEY_MINIMAP2_OPTIONS]
    log: os.path.join(config[KEY_TEMPDIR],"logs","{barcode}.minimap2_initial.log")
    output:
        paf = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","filtered_reads.paf")
    shell:
        """
        minimap2 -t {threads} {params.minimap2_options} --secondary=no --paf-no-hit \
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
        hits = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","hits_reads.csv"),
        multi_out = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","multi_mapped_refs.csv"),
        parsing = os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","hits_filter_description.csv")
    run:
        parse_paf_file(input.map_file,
                        output.csv,
                        output.parsing,
                        output.hits,
                        output.multi_out,
                        config[KEY_REFERENCE_SEQUENCES],
                        config[KEY_POSITIVE_REFERENCES],
                        config[KEY_INCLUDE_POSITIVE_REFERENCES],
                        params.barcode,
                        config[KEY_ANALYSIS_MODE],
                        config[KEY_MIN_MAP_QUALITY],
                        config[KEY_REFERENCE_GROUP_FIELD],
                        config)

rule write_hit_fastq:
    input:
        csv = rules.assess_broad_diversity.output.csv,
        hits = rules.assess_broad_diversity.output.hits,
        fastq = rules.filter_by_length.output.fastq
    params:
        primer_len = config[KEY_PRIMER_LENGTH],
        barcode = "{barcode}",
        outdir = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups"),
        publish_dir = os.path.join(config[KEY_OUTDIR],"published_data","{barcode}")
    output:
        txt = os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt"),
        pub_txt = os.path.join(config[KEY_OUTDIR],"published_data","{barcode}","prompt.txt")
    run:
        shell("touch {output.txt:q} && touch {output.pub_txt:q}")
        print(green(params.barcode))
        to_write = write_out_fastqs(input.csv,input.hits,input.fastq,params.outdir,params.primer_len,config)

        for ref in to_write:
            written = os.path.join(params.outdir,f"{ref}.fastq")
            published = os.path.join(params.publish_dir,f"{ref}.fastq")
            shell(f"cp '{written}' '{published}'")

        write_out_ref_fasta(to_write,config[KEY_REFERENCE_SEQUENCES],params.outdir)


rule gather_diversity_report:
    input:
        refs = expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","initial_processing","refs_present.csv"), barcode=config[KEY_BARCODES]),
        txt = expand(os.path.join(config[KEY_TEMPDIR],"{barcode}","reference_groups","prompt.txt"), barcode=config[KEY_BARCODES]),
        ref = config[KEY_REFERENCE_SEQUENCES]
    output:
        refs= os.path.join(config[KEY_TEMPDIR],SAMPLE_COMPOSITION),
        summary = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY),
        yaml = os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG)
    run:
        barcode_config = diversity_report(input.refs,output.refs,output.summary,input.ref,config)
        
        with open(output.yaml, 'w') as fw:
            yaml.dump(barcode_config, fw) 
        print(yellow("-----------------------"))
