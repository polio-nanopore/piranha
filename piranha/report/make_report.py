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

def make_sample_report(report_to_generate,
                        variation_file,
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
    
    with open(variation_file,"r") as f:
        variation_json = json.load(f)
        for reference in variation_json:

            all_sites_data = variation_json[reference]["variation"]
            
            coocc_data = variation_json[reference]["coocc"]
            if coocc_data and "Sabin" in reference:
                data_for_report[reference][KEY_COOCCURRENCE_INFO] = coocc_data

            snp_sites = data_for_report[reference][KEY_SNP_SITES]
            masked_sites = data_for_report[reference][KEY_MASKED_SITES]
            indel_sites = data_for_report[reference][KEY_INDEL_SITES]

            

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

def make_detailed_csv(data_for_report,barcodes_csv,output):
    """
    SUMMARY TABLE
        list of: {KEY_BARCODE:record_barcode,
                    KEY_SAMPLE:record_sample,
                    "Sample call": call,
                    KEY_REFERENCE_GROUP:reference_group,
                    "Number of mutations": int(var_count),
                    }
    COMPOSITION TABLE
        rows of:
        sample,barcode,Sabin1-related,Sabin2-related,Sabin3-related,WPV1,WPV2,WPV3,NonPolioEV,unmapped
        MixedTest,barcode01,0,246,250,0,0,0,1,0
        PureTest,barcode02,0,709,0,0,0,0,7,12
        WTTest,barcode03,0,1,3,246,0,0,1,10
        VDPVTest,barcode04,0,0,250,0,0,0,0,0
        negative,barcode05,0,0,5,0,0,0,0,0
        positively,barcode06,0,0,5,0,0,0,0,0
        somemixed,barcode07,0,138,339,9,0,0,74,23

    """

    #detailed report cols: "sample,barcode,EPID,Institute,...[from_csv]...,[detailed composition header fields],..."

    """
    #detailed header fields
    "Sabin-related|closest_reference","Sabin-related|num_reads","Sabin-related|nt_diff_from_reference","Sabin-related|pcent_match","Sabin-related|classification",
    "WT-Polio|closest_reference","WT-Polio|num_reads","WT-Polio|nt_diff_from_reference","WT-Polio|pcent_match","WT-Polio|classification",
    "NonPolioEV|num_reads","comments"
    """

    combined_barcodes = collections.defaultdict(dict)
    for row in data_for_report[KEY_SUMMARY_TABLE]:
        combined_barcodes[row[KEY_BARCODE]][row[KEY_REFERENCE_GROUP]] = row

    sample_metadata = {}
    metadata_headers = []
    with open(barcodes_csv,"r") as f:
        reader = csv.DictReader(f)
        metadata_headers = reader.fieldnames
        for row in reader:
            sample_metadata[row[KEY_BARCODE]] = row

    header_fields = [KEY_SAMPLE,KEY_BARCODE,KEY_EPID,KEY_INSTITUTE]
    for i in metadata_headers:
        if i not in header_fields:
            header_fields.append(i)

    for i in DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS:
        header_fields.append(i)

    with open(output,"w") as fw:
        writer = csv.DictWriter(fw,fieldnames=header_fields,lineterminator="\n")
        writer.writeheader()
        

        for row in data_for_report[KEY_COMPOSITION_TABLE]:
            to_write = {}
            to_write = {}
            for i in header_fields:
                to_write[i] = ""
            to_write[KEY_SAMPLE] = row[KEY_SAMPLE]
            to_write[KEY_BARCODE] = row[KEY_BARCODE]
            
            metadata = sample_metadata[row[KEY_BARCODE]]
            
            #take in info from barcodes csv here
            for field in metadata:
                to_write[field] = metadata[field]

            # read count
            for col in row:
                if col not in [KEY_SAMPLE,KEY_BARCODE,"unmapped"]:
                    to_write[f"{col}|num_reads"] = row[col]

            for ref_group in combined_barcodes[row[KEY_BARCODE]]:
                info = combined_barcodes[row[KEY_BARCODE]][ref_group]
                to_write[f"{ref_group}|closest_reference"] = info[KEY_REFERENCE]
                to_write[f"{ref_group}|nt_diff_from_reference"] = info["Number of mutations"]
                to_write[f"{ref_group}|pcent_match"] = info[KEY_PERCENT]
                to_write[f"{ref_group}|classification"] = info["Sample classification"]
        
            writer.writerow(to_write)

def well_to_dict(barcode_well_map,j,i,c):
    if i<10:
        well = f"{j}0{i}"
        
    else:
        well = f"{j}{i}"
        
    if c<10:
        barcode_well_map[well] = f"barcode0{c}"
    else:
        barcode_well_map[well] = f"barcode{c}"
        
def assign_bcode_to_well(orientation):
    barcode_well_map = {}
    c = 0
    if orientation == "vertical":
        for i in range(1,13):
            for j in ["A","B","C","D","E","F","G","H"]:
                c +=1
                well_to_dict(barcode_well_map,j,i,c)
                
    elif orientation == "horizontal":
        for j in ["A","B","C","D","E","F","G","H"]:
            for i in range(1,13):
                c +=1
                well_to_dict(barcode_well_map,j,i,c)
                
                
    return barcode_well_map
                

def barcode_to_well(barcode_csv,orientation):
    
    barcode_well_map = {}
    with open(barcode_csv,"r") as f:
        reader = csv.DictReader(f)
        if "well" in reader.fieldnames:
            for row in reader:
                barcode_well_map[row["well"]]=row["barcode"]
        else:
            barcode_well_map = assign_bcode_to_well(orientation)
            
        return barcode_well_map
            
        


def data_for_plate_viz(positives_for_plate_viz,barcode_csv,orientation):
    
    
    barcode_well_map = barcode_to_well(barcode_csv,orientation)   
    wells_to_json = []
    
    
    all_positive_types = set()
    all_positive_types.add("All")
    for i in positives_for_plate_viz:
        types = positives_for_plate_viz[i].keys()
        for t in types:
            all_positive_types.add(t)
    all_positive_types = sorted(list(all_positive_types))
            
    
    for j in ["A","B","C","D","E","F","G","H"]:
        for i in range(1,13):
            if i<10:
                well = f"{j}0{i}"
            else:
                well = f"{j}{i}"
            info = {
                "x":i,
                "y":j,
                "All":"Absent"
            }
            
            pos_type = ""

            if well in barcode_well_map:
                
                barcode = barcode_well_map[well]
                info["Barcode"] = barcode
                
                
                if barcode in positives_for_plate_viz:
                    positives = positives_for_plate_viz[barcode]
                    
                    for i in all_positive_types:
                        if i in positives:
                            info["All"] = "Present"
                            info[i] = "Present"
                        else:
                            info[i] = "Absent"
                else:
                    for i in all_positive_types:
                        info[i] = "Absent"
            else:
                info["Barcode"] = ""
                
            wells_to_json.append(info)
            
    return json.dumps(wells_to_json), all_positive_types

def make_output_report(report_to_generate,barcodes_csv,preprocessing_summary,sample_composition,consensus_seqs,detailed_csv_out,config):
    
    # which are the negative controls and positive controls
    negative_control = config[KEY_NEGATIVE]
    positive_control = config[KEY_POSITIVE]

    # which will be flagged as high npev samples
    flagged_high_npev = []

    # does it pass control or not
    control_status = {negative_control:True,positive_control:True}
    #are there any controls
    show_control_table = False
    
    # collate data for tables in the report
    data_for_report = {KEY_SUMMARY_TABLE:[],KEY_COMPOSITION_TABLE:[]}
    positives_for_plate_viz = collections.defaultdict(dict)

    with open(preprocessing_summary,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # add directly to data for report composition table (2)
            data_for_report[KEY_COMPOSITION_TABLE].append(row)

            # handling controls
            if row[KEY_SAMPLE] == negative_control:
                show_control_table = True
                for col in row:
                    if col not in [KEY_BARCODE,KEY_SAMPLE]:
                        # if there are any reads above the min read threshold
                        # sample fails as a negative control
                        if int(row[col])>config[KEY_MIN_READS]:
                            control_status[negative_control] = False
            elif row[KEY_SAMPLE] == positive_control:
                show_control_table = True
                # if the npev reads are below the min read threshold
                # sample fails as a positive control
                if int(row["NonPolioEV"])<config[KEY_MIN_READS]:
                    control_status[positive_control] = False
            else:
                # deciding whether to flag the npev table
                # if there are any polio reads
                # if pcent of sample that is npev greater than min pcent
                # flag the table

                total_reads = 0
                total_polio_reads = 0
                for col in row:
                    if col not in [KEY_SAMPLE,KEY_BARCODE]:
                        total_reads+=int(row[col])

                        if int(row[col])>config[KEY_MIN_READS]:
                            positives_for_plate_viz[row[KEY_BARCODE]][col] = int(row[col])
                            
                        if col not in ["NonPolioEV","unmapped"]:
                            total_polio_reads+= int(row[col])

                if total_polio_reads>0:
                    proportion_npev = 100*(int(row["NonPolioEV"])/total_reads)
                    if proportion_npev > config[KEY_MIN_PCENT]:
                        flagged_high_npev.append(row[KEY_SAMPLE])

    # to check if there are identical seqs in the run
    identical_seq_check = collections.defaultdict(list)

    plate_json, positive_types = data_for_plate_viz(positives_for_plate_viz,barcodes_csv,config[KEY_ORIENTATION])

    for record in SeqIO.parse(consensus_seqs,KEY_FASTA):
        identical_seq_check[str(record.seq)].append(record.description)
        
        fields = record.description.split("|")

        # the first fields there will always be
        record_sample,record_barcode,reference_group,reference,var_count,var_string = fields[:6]
        # additional fields that have been added to header
        additional_fields = fields[6:]

        additional_info = {}
        for field in additional_fields:
            key,value = field.split("=")
            additional_info[key]=value
        
        length_of_seq = len(record)

        prop_diff = int(var_count)/length_of_seq
        pcent_match = round((1-prop_diff)*100, 2)

        call = reference_group
        if reference_group.startswith("Sabin"):
            # configured number of mutations in sabin for the call threshold of VDPV
            call_threshold = CALL_THRESHOLD_DICT[reference_group]
            if int(var_count) > call_threshold:
                call = "VDPV"
            elif var_count == 0:
                call = "Sabin"
            else:
                call = "Sabin-like"

            info = {KEY_BARCODE:record_barcode,
                KEY_SAMPLE:record_sample,
                "Sample classification": call,
                KEY_REFERENCE_GROUP:reference_group,
                KEY_REFERENCE:reference,
                "Number of mutations": int(var_count),
                KEY_PERCENT:pcent_match
                }
        else:
            info = {KEY_BARCODE:record_barcode,
                KEY_SAMPLE:record_sample,
                "Sample classification": call,
                KEY_REFERENCE_GROUP:reference_group,
                KEY_REFERENCE:reference,
                "Number of mutations": "NA",
                KEY_PERCENT:pcent_match
                }

        data_for_report[KEY_SUMMARY_TABLE].append(info)

    # flag if there are identical seqs
    flagged_seqs = []
    for seq in identical_seq_check:
        if len(identical_seq_check[seq]) > 1:
            flagged_seqs.append(identical_seq_check[seq])

    # do we need the control table
    data_for_report[KEY_CONTROL_STATUS] = control_status
    
    # composition table header
    if config[KEY_ANALYSIS_MODE] == VALUE_ANALYSIS_MODE_WG:
        config[KEY_COMPOSITION_TABLE_HEADER] = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG
    else:
        config[KEY_COMPOSITION_TABLE_HEADER] = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1

    # summary table header
    config[KEY_SUMMARY_TABLE_HEADER] = SAMPLE_SUMMARY_TABLE_HEADER_FIELDS
    
    # detailed csv for download (1)
    make_detailed_csv(data_for_report,barcodes_csv,detailed_csv_out)

    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(),
                    run_name=config[KEY_RUN_NAME],
                    version = __version__,
                    show_control_table = show_control_table,
                    plate_json = plate_json,
                    positive_types = positive_types,
                    data_for_report = data_for_report,
                    flagged_seqs = flagged_seqs,
                    detailed_csv_out = detailed_csv_out,
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
