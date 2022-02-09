import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.analysis.parse_paf import parse_paf_file
from piranha.analysis.filter_lengths import filter_reads_by_length
from piranha.report.make_report import make_report

"""
requires snakemake --snakefile piranha/scripts/report_runner.smk --cores 1 --config barcode=barcode01 outdir=analysis_2022-02-09
"""
barcode = config['barcode']

rule generate_report:
    input:
        haplotype_info = os.path.join(config[KEY_OUTDIR],f"{barcode}","processed_data","haplotypes.csv"),
        read_lengths = os.path.join(config[KEY_OUTDIR],f"{barcode}","processed_data","lengths.txt"),
        variation_info = os.path.join(config[KEY_OUTDIR],f"{barcode}","processed_data","variation_info.json"),
        yaml = os.path.join(config[KEY_OUTDIR],f"{barcode}","processed_data","haplotype_config.yaml")
    output:
        html = os.path.join(config[KEY_OUTDIR],f"{barcode}","analysis_report.html")
    run:
        with open(input.yaml, 'r') as f:
            config_loaded = yaml.safe_load(f)

        make_report(output.html,input.read_lengths,input.variation_info,input.haplotype_info,barcode,config_loaded)
