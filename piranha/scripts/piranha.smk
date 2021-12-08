import os
import collections
from Bio import SeqIO
import yaml

import apollofunks as qcfunk

if config.get("force"):
    config["force"] = "--forceall "

##### Target rules #####

rule all:
    input:
        os.path.join(config["outdir"], "reports","report.csv")

rule demultiplex:
    params:
        outdir = os.path.join(config["outdir"],"demultiplexed_reads")
    threads:
        workflow.cores
    output:
        demux_prompt = os.path.join(config["outdir"],"demultiplexed_reads", "demuxed.txt")
    run:
        if config["demultiplex"]:
            shell("""
            {config[path_to_guppy]:q} \
            -i {config[read_path]:q} \
            -s {params.outdir:q} \
            -t {threads} \
            --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" \
            && touch {output.demux_prompt:q}
            """)
            config["read_path"] = params.outdir
        else:
            shell("touch {output.demux_prompt:q}")

rule count_cpgs:
    input:
        snakefile = os.path.join(workflow.current_basedir,"count_cpgs.smk"),
        demux_prompt = os.path.join(config["outdir"],"demultiplexed_reads", "demuxed.txt"),
        genes = config["genes"],
        cpg_sites = config["cpg_sites"],
        primer_sequences = config["primer_sequences"],
        matrix_file = config["matrix_file"]
    threads:
        workflow.cores
    output:
        yaml = os.path.join(config["outdir"], "config.yaml"),
        counts = os.path.join(config["outdir"], "reports","cpg_counts.csv"),
        wide = os.path.join(config["outdir"], "reports","cpg_wide.csv")
    run:
        print(qcfunk.green("Barcodes found:"))

        if config["demultiplex"]:
            config["read_path"] = os.path.join(config["outdir"],"demultiplexed_reads")

        barcodes = []
        for r,d,f in os.walk(config["read_path"]):
            for name in d:
                if name.startswith("barcode"):
                    print(name)
                    barcodes.append(name)

        barcode_str = ','.join(barcodes)
        config["barcodes"] = barcodes
        
        with open(output.yaml, 'w') as fw:
            yaml.dump(config, fw) #so at the moment, every config option gets passed to the report

        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                                "{config[force]} "
                                "{config[log_string]} "
                                "--directory {config[tempdir]:q} "
                                "--configfile {output.yaml:q} " 
                                "--cores {threads} ")
