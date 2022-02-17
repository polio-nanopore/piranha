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

def make_sample_report(report_to_generate,variation_file,consensus_seqs,barcode,config):

    references = config[barcode]

    data_for_report = {}

    for reference in references:
        data_for_report[reference] = {}
        data_for_report[reference]["snipit_svg"] = get_snipit(reference,os.path.join(config[KEY_TEMPDIR],"snipit",f"{reference}.svg"))

    info = {}

    for record in SeqIO.parse(consensus_seqs,"fasta"):
        #f"{sample]}|{barcode}|{reference_group}|{var_count}|{var_string}|{date}"
        record_sample,record_barcode,reference_group,ref,var_count,var_string,collection_date = record.id.split("|")
        if ref == reference:
            info = {KEY_BARCODE:barcode,
                    KEY_SAMPLE:record_sample,
                    KEY_REFERENCE_GROUP:reference_group,
                    "Number of mutations": var_count,
                    "Variants":var_string,
                    "Collection date":date
                    }

    data_for_report[reference]["summary_data"] = info

    with open(variation_file,"r") as f:
        var_data = json.load(f)
        for reference in var_data:
            data_for_report[reference]["variation_info"] = var_data[reference]
            
    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(), 
                    version = __version__,
                    barcode = barcode,
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

    data_for_report = {"summary_table":[]}
    with open(preprocessing_summary,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report["summary_table"].append(row)
    
    config["table_header"] = SAMPLE_SUMMARY_HEADER_FIELDS
    
    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(),
                    run_name=config[KEY_RUN_NAME],
                    version = __version__,
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
