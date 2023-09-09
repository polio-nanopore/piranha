import os
import collections
from Bio import SeqIO
import yaml
import pkg_resources
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
from piranha.report.make_report import make_output_report

"""
requires snakemake --snakefile piranha/scripts/report_runner.smk --cores 1 --config outdir=analysis_2022-02-09
"""
package_datafile = os.path.join("data","report.mako")
data = pkg_resources.resource_filename('piranha',package_datafile)
config[KEY_REPORT_TEMPLATE] = data

rule generate_report:
    input:
        summary_csv=os.path.join(config[KEY_OUTDIR],PREPROCESSING_SUMMARY),
        composition_csv=os.path.join(config[KEY_OUTDIR],SAMPLE_COMPOSITION),
        yaml = os.path.join(config[KEY_OUTDIR],PREPROCESSING_CONFIG),
        seqs = os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS),
        background_data = os.path.join(config[KEY_OUTDIR],"phylogenetics","annotations.csv"),
        detailed_csv = os.path.join(config[KEY_OUTDIR],"detailed_run_report.csv")
    output:
        report =os.path.join(config[KEY_OUTDIR],OUTPUT_REPORT)
    run:
        with open(input.yaml,"r") as f:
            config_loaded = yaml.safe_load(f) 
        
        make_output_report(output.report,config_loaded["barcodes_csv"],input.summary_csv,input.composition_csv,input.seqs,input.detailed_csv,input.background_data,config_loaded)