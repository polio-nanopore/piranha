import json
import csv
import collections
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date

from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import os
from Bio import SeqIO

from piranha import __version__
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


def get_snipit(reference,snipit_file):
    snipit_svg = ""
    with open(snipit_file,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            snipit_svg+=f"{l}\n"
    return snipit_svg

def make_sample_report(report_to_generate,variation_file,consensus_seqs,masked_variants,barcode,config):

    references = config[barcode]
    
    data_for_report = {}

    for reference in references:
        data_for_report[reference] = {}
        data_for_report[reference]["snipit_svg"] = get_snipit(reference,os.path.join(config[KEY_TEMPDIR],f"{barcode}","snipit",f"{reference}.svg"))

    info_dict = {}

    for record in SeqIO.parse(consensus_seqs,"fasta"):
        #f"{sample]}|{barcode}|{reference_group}|{var_count}|{var_string}|{date}"
        record_sample,record_barcode,reference_group,reference,var_count,var_string,collection_date = record.description.split("|")
        if barcode == record_barcode:
            sample = record_sample
            info = {KEY_BARCODE:barcode,
                    KEY_SAMPLE:record_sample,
                    KEY_REFERENCE_GROUP:reference_group,
                    "Number of mutations": int(var_count),
                    "Variants":var_string,
                    "Collection date":date
                    }
            info_dict[reference] = info

            data_for_report[reference]["snp_sites"] = []
            data_for_report[reference]["indel_sites"] = []

            for var in var_string.split(";"):
                site = var.split(":")[0]
                if "ins" in var or "del" in var:
                    data_for_report[reference]["indel_sites"].append(int(site))
                else:
                    data_for_report[reference]["snp_sites"].append(int(site))

    for reference in info_dict:
        data_for_report[reference]["masked_sites"] = []
        data_for_report[reference]["summary_data"] = info_dict[reference]

    with open(masked_variants,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report[row["reference"]]["masked_sites"].append(int(row["site"]))

    with open(variation_file,"r") as f:
        var_data = json.load(f)
        for reference in var_data:
            data_for_report[reference]["variation_info"] = var_data[reference]
    
    template_dir = os.path.abspath(os.path.dirname(config[KEY_BARCODE_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_BARCODE_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(), 
                    version = __version__,
                    barcode = barcode,
                    sample = sample,
                    data_for_report = data_for_report,
                    config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print(green("Generating: ") + f"{report_to_generate}")
        fw.write(buf.getvalue())

def make_output_report(report_to_generate,preprocessing_summary,sample_composition,config):

    control_status = {"negative":True,"positive":True}
    data_for_report = {"summary_table":[]}
    show_control_table = False
    with open(preprocessing_summary,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report["summary_table"].append(row)
            if row["sample"] == "negative":
                show_control_table = True
                for col in row:
                    if col not in ["barcode","sample"]:
                        if int(row[col])>config[KEY_MIN_READS]:
                            control_status["negative"] = False
            elif row["sample"] == "positive":
                show_control_table = True
                if int(row["NonPolioEV"])<config[KEY_MIN_READS]:
                    control_status["positive"] = False
                    

    data_for_report["control_status"] = control_status

    config["table_header"] = SAMPLE_SUMMARY_HEADER_FIELDS
    
    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(),
                    run_name=config[KEY_RUN_NAME],
                    version = __version__,
                    show_control_table = show_control_table,
                    data_for_report = data_for_report,
                    config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print(green("Generating: ") + f"{report_to_generate}")
        fw.write(buf.getvalue())
