import os
import collections
from Bio import SeqIO
import yaml
import importlib.resources
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.report.make_report import make_output_report

"""
requires snakemake --snakefile piranha/scripts/report_runner.smk --cores 1 --config outdir=analysis_2023-11-14
"""
package_datafile = os.path.join("data","report.mako")
    resource = importlib.resources.files("piranha").joinpath("data", package_datafile)
    with importlib.resources.as_file(resource) as data_path:
        data = str(data_path)
config[KEY_REPORT_TEMPLATE] = data

rule generate_report:
    input:
        summary_csv=os.path.join(config[KEY_OUTDIR],PREPROCESSING_SUMMARY),
        composition_csv=os.path.join(config[KEY_OUTDIR],SAMPLE_COMPOSITION),
        yaml = os.path.join(config[KEY_OUTDIR],PREPROCESSING_CONFIG),
        cns_yaml = os.path.join(config[KEY_OUTDIR],"consensus_config.yaml"),
        seqs = os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS),
        detailed_csv = os.path.join(config[KEY_OUTDIR],"detailed_run_report.csv")
    output:
        report =os.path.join(config[KEY_OUTDIR],OUTPUT_REPORT)
    run:
        with open(input.yaml,"r") as f:
            config_loaded = yaml.safe_load(f) 
        with open(input.cns_yaml, 'r') as f:
            cns_config_loaded = yaml.safe_load(f)
        make_output_report(output.report,config_loaded["barcodes_csv"],input.summary_csv,input.composition_csv,input.seqs,input.detailed_csv,input.background_data,cns_config_loaded,config_loaded)