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


def get_snipit(taxon,snipit_file):
    snipit_svg = ""
    with open(snipit_file,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            snipit_svg+=f"{l}\n"
    return snipit_svg

def make_report(report_to_generate,read_length_file,variation_file,haplotype_info,barcode,config):

    taxa_present = config[KEY_TAXA_PRESENT]

    data_for_report = {}
    data_for_report["lengths"] = []
    with open(read_length_file, "r") as f:
        for l in f:
            l = l.rstrip("\n")
            try:
                data_for_report["lengths"].append({"x":int(l)})
            except:
                pass

    for taxon in taxa_present:
        data_for_report[taxon] = {}
        data_for_report[taxon]["snipit_svg"] = get_snipit(taxon,os.path.join(config[KEY_TEMPDIR],"snipit",f"{taxon}.svg"))

    summary_data = []
    num_sites = {}
    with open(haplotype_info,"r") as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            new_row = row
            num_sites[new_row["taxon"]] = len(new_row["sites"].split(";"))
            summary_data.append(new_row)

    data_for_report["summary_data"] = summary_data

    with open(variation_file,"r") as f:
        var_data = json.load(f)
        for taxon in var_data:
            data_for_report[taxon]["variation_info"] = var_data[taxon]
            data_for_report[taxon]["num_sites"] = num_sites[taxon]
            
    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(), 
                    version = __version__,
                    barcode = barcode,
                    data_for_report = data_for_report,
                    taxa_present = config["taxa_present"],
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