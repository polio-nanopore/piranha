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

# def include_call_info():

def make_sample_report(report_to_generate,
                        variation_file,
                        co_occurrence_file,
                        consensus_seqs,
                        masked_variants,
                        barcode,
                        config):

    references = config[barcode]
    
    data_for_report = {}

    for reference in references:
        data_for_report[reference] = {}
        data_for_report[reference][KEY_SNIPIT_SVG] = get_snipit(reference,os.path.join(config[KEY_TEMPDIR],f"{barcode}","snipit",f"{reference}.svg"))

    info_dict = {}
    sequences = ""


    for record in SeqIO.parse(consensus_seqs,KEY_FASTA):
            
        
        fields = record.description.split("|")
        
        record_sample,record_barcode,reference_group,reference,var_count,var_string = fields[:6]
        additional_fields = fields[6:]

        
        additional_info = {}
        for field in additional_fields:
            key,value = field.split("=")
            additional_info[key]=value

        if barcode == record_barcode:
            
            sequences+= f">{record.description}<br>{record.seq}<br>"
            
            sample = record_sample
            info = {KEY_BARCODE:barcode,
                    KEY_SAMPLE:record_sample,
                    KEY_REFERENCE_GROUP:reference_group,
                    "Number of mutations": int(var_count),
                    "Variants":var_string
                    }

            for key in additional_info:
                info[key] = additional_info[key]

            info_dict[reference] = info

            data_for_report[reference][KEY_SNP_SITES] = []
            data_for_report[reference][KEY_INDEL_SITES] = []

            for var in var_string.split(";"):
                site = var.split(":")[0]
                if "ins" in var or "del" in var:
                    data_for_report[reference][KEY_INDEL_SITES].append(int(site))
                else:
                    try:
                        site = int(site)
                        data_for_report[reference][KEY_SNP_SITES].append(site)
                    except:
                        data_for_report[reference][KEY_SNP_SITES].append(site)

    

    for reference in info_dict:
        data_for_report[reference][KEY_MASKED_SITES] = []
        data_for_report[reference][KEY_SUMMARY_DATA] = info_dict[reference]

    with open(masked_variants,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report[row[KEY_REFERENCE]][KEY_MASKED_SITES].append(int(row[KEY_SITE]))
    
    with open(co_occurrence_file,"r") as f:
        co_occurrence_data = json.load(f)
        for reference in co_occurrence_data:
            comb_data = []
            #{"Poliovirus3-Sabin_AY184221": {"17T;160A": 42, "17T;160G": 48}}
            for comb in co_occurrence_data[reference]:
                comb_info = {}
                comb_info["Combination"] = comb
                comb_info["Percentage"] = co_occurrence_data[reference][comb]
                comb_data.append(comb_info)

            data_for_report[reference][KEY_COOCCURRENCE_INFO] = comb_data

    with open(variation_file,"r") as f:
        var_data = json.load(f)
        for reference in var_data:
            snp_sites = data_for_report[reference][KEY_SNP_SITES]
            masked_sites = data_for_report[reference][KEY_MASKED_SITES]
            indel_sites = data_for_report[reference][KEY_INDEL_SITES]

            all_sites_data = var_data[reference]

            annotated_site_data = []
            for site_data in all_sites_data:
                if site_data["Position"] in snp_sites:
                    site_data["snp_type"] = "SNP"
                    site_data["colour"] = "#133239"
                    site_data["size"] = 20
                elif site_data["Position"] in masked_sites:
                    site_data["snp_type"] = "Masked Variant"
                    site_data["colour"] = "#e68781"
                    site_data["size"] = 20
                elif site_data["Position"] in indel_sites:
                    site_data["snp_type"] = "Insertion/ Deletion"
                    site_data["colour"] = "#B99C0C"
                    site_data["size"] = 20
                else:
                    site_data["snp_type"] = "Background variation"
                    site_data["colour"] = "#ACAFB0"
                    site_data["size"] = 10

                annotated_site_data.append(site_data)

            data_for_report[reference][KEY_VARIATION_INFO] = annotated_site_data

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
                    sequences = sequences,
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

def make_output_report(report_to_generate,preprocessing_summary,sample_composition,consensus_seqs,config):
    negative_control = config[KEY_NEGATIVE]
    positive_control = config[KEY_POSITIVE]

    flagged_high_npev = []

    control_status = {negative_control:True,positive_control:True}
    data_for_report = {KEY_SUMMARY_TABLE:[],KEY_COMPOSITION_TABLE:[]}
    show_control_table = False
    with open(preprocessing_summary,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data_for_report[KEY_COMPOSITION_TABLE].append(row)
            if row[KEY_SAMPLE] == negative_control:
                show_control_table = True
                for col in row:
                    if col not in [KEY_BARCODE,KEY_SAMPLE]:
                        if int(row[col])>config[KEY_MIN_READS]:
                            control_status[negative_control] = False
            elif row[KEY_SAMPLE] == positive_control:
                show_control_table = True
                if int(row["NonPolioEV"])<config[KEY_MIN_READS]:
                    control_status[positive_control] = False
            else:
                if int(row["NonPolioEV"])>100:
                    flagged_high_npev.append(row[KEY_SAMPLE])

    identical_seq_check = collections.defaultdict(list)

    for record in SeqIO.parse(consensus_seqs,KEY_FASTA):
        identical_seq_check[str(record.seq)].append(record.description)

        fields = record.description.split("|")

        record_sample,record_barcode,reference_group,reference,var_count,var_string = fields[:6]
        additional_fields = fields[6:]

        additional_info = {}
        for field in additional_fields:
            key,value = field.split("=")
            additional_info[key]=value


        call = reference_group
        if reference_group.startswith("Sabin"):
            call_threshold = CALL_THRESHOLD_DICT[reference_group]
            if int(var_count) > call_threshold:
                call = "VDPV"
            elif var_count == 0:
                call = "Sabin"
            else:
                call = "Sabin-like"

            info = {KEY_BARCODE:record_barcode,
                KEY_SAMPLE:record_sample,
                "Sample call": call,
                KEY_REFERENCE_GROUP:reference_group,
                "Number of mutations": int(var_count),
                }
        else:
            info = {KEY_BARCODE:record_barcode,
                KEY_SAMPLE:record_sample,
                "Sample call": call,
                KEY_REFERENCE_GROUP:reference_group,
                "Number of mutations": "NA",
                }


        
                
        data_for_report[KEY_SUMMARY_TABLE].append(info)

    flagged_seqs = []
    for seq in identical_seq_check:
        if len(identical_seq_check[seq]) > 1:
            flagged_seqs.append(identical_seq_check[seq])

    data_for_report[KEY_CONTROL_STATUS] = control_status
    
    if config[KEY_ANALYSIS_MODE] == VALUE_ANALYSIS_MODE_WG_2TILE:
        config[KEY_COMPOSITION_TABLE_HEADER] = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG
    else:
        config[KEY_COMPOSITION_TABLE_HEADER] = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1

    config[KEY_SUMMARY_TABLE_HEADER] = SAMPLE_SUMMARY_TABLE_HEADER_FIELDS
    
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
                    flagged_seqs = flagged_seqs,
                    flagged_high_npev = flagged_high_npev,
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
