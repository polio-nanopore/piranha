#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv
import json
from itertools import groupby

from piranha.utils.config import *

def gather_fasta_files(summary_info, barcodes_csv, input_cns_list,all_metdata,output_file,publish_dir,config):
    if not os.path.exists(publish_dir):
        os.mkdir(publish_dir)
    
    analysis_info = collections.defaultdict(list)
    with open(summary_info, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            analysis_info[row[KEY_BARCODE]].append(row)

    input_metadata = {}
    handle_dict = {}
    with open(barcodes_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            barcode = row[KEY_BARCODE]
            input_metadata[barcode] = row
            
            if barcode not in handle_dict:
                if not os.path.exists(os.path.join(publish_dir, f"{barcode}")):
                    os.mkdir(os.path.join(publish_dir, f"{barcode}"))
                handle_dict[barcode] = open(os.path.join(publish_dir, f"{barcode}",f"{barcode}.consensus.fasta"),"w")

    fasta_header_fields = config[KEY_FASTA_HEADER_FIELDS]

    with open(output_file,"w") as fw:
        for cns_file in input_cns_list:
            for record in SeqIO.parse(cns_file, KEY_FASTA):
                if record:
                    cns_info= record.description.split(" ")
                    ref_hap,barcode,var_count,var_string=cns_info[0].split("|")
                    
                    # for parsing haplotypes
                    ref_list = ref_hap.split(".")
                    ref = ".".join(ref_list[:-1])
                    hap = ref_list[-1]

                    metadata = input_metadata[barcode]
                    metadata[KEY_VARIANT_COUNT] = var_count
                    metadata[KEY_VARIANTS] = var_string
                    metadata[KEY_CNS_ID] = hap
                    metadata[KEY_RUNNAME] = config[KEY_RUNNAME]
                    metadata[KEY_RUNID] = config[KEY_RUNID]

                    info = []
                    for row in analysis_info[barcode]:
                        if row[KEY_REFERENCE] == ref:
                            info = row

                    for col in info:
                        metadata[col] = info[col]

                    record_id = f""
                    for field in config[KEY_FASTA_HEADER_FIELDS]:
                        
                        if field in metadata:
                            record_id += f"{metadata[field]}|"
                            print(field, metadata[field])
                        else:
                            record_id += f"|"
                            print(field, "missing")

                    # record_id += f"|{hap}|{metadata[KEY_REFERENCE_GROUP]}"

                    # if KEY_EPID in metadata:
                    #     record_id += f"|{metadata[KEY_EPID]}"
                    # else:
                    #     record_id += "|"

                    # if KEY_DATE in metadata:
                    #     record_id += f"|{metadata[KEY_DATE]}"
                    # else:
                    #     record_id += "|"

                    record_id += f" {KEY_BARCODE}={barcode}"
                    record_id += f" {KEY_REFERENCE}={ref}"
                    record_id += f" {VALUE_REFERENCE_GROUP_FIELD}={info[KEY_REFERENCE_GROUP]}"

                    # if runname:
                    #     record_id += f" {KEY_RUNNAME}={runname}"

                    if "Sabin" in ref:
                        record_id += f" {KEY_VARIANT_COUNT}={var_count}"
                        record_id += f" {KEY_VARIANTS}={var_string}"
                    else:
                        record_id += f" {KEY_VARIANT_COUNT}=NA"
                        record_id += f" {KEY_VARIANTS}=NA"

                    if all_metdata:
                        
                        for col in metadata:
                            if col != KEY_SAMPLE and col != KEY_BARCODE:
                                record_id += f" {col}={metadata[col]}"
                    """
                    record header is:
                    >SAMPLE|REFERENCE_GROUP|CNS_ID|EPID|DATE barcode=barcode01 variant_count=8 variants=17:CT;161:CT;427:GA;497:AC;507:CT;772:AG;822:CT;870:CA 

                    if "all_metadata" then everything else gets added to the description
                    """
                    fw.write(f">{record_id}\n{record.seq}\n")
                    handle_dict[barcode].write(f">{record_id}\n{record.seq}\n")
    
    for handle in handle_dict:
        handle_dict[handle].close()

def get_sample(barcodes_csv,barcode):
    sample = ""
    with open(barcodes_csv,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row[KEY_BARCODE] == barcode:
                    sample = row[KEY_SAMPLE]
    return sample




def group_consecutive_sites(lst):
    out = []
    for _, g in groupby(enumerate(lst), lambda k: k[0] - k[1]):
        start = next(g)[1]
        end = list(v for _, v in g) or [start]
        out.append(range(start, end[-1] + 1))
    return out

def get_mask_dict(mask_file):
    mask_dict = {}
    with open(mask_file,"r") as f:
        mask_json = json.load(f)
        
        for ref in mask_json:

            sites_to_mask = sorted(mask_json[ref])
            mask_ranges = group_consecutive_sites(sites_to_mask)
            mask_dict[ref] = mask_ranges
    return mask_dict
            

def mask_low_coverage(mask_file, sequences,output):
    mask_dict = get_mask_dict(mask_file)
    records = 0
    with open(output,"w") as fw:
        for record in SeqIO.parse(sequences,"fasta"):
            seq_id = record.id.split("|")[0]
            if seq_id in mask_dict and mask_dict[seq_id]:
                mask_ranges = mask_dict[seq_id]
                n_count = 0
                new_seq = str(record.seq)
                for site in mask_ranges:

                    start,stop = site[0]-1,site[-1]+1

                    if site[0] == 0:
                        length = stop
                        n_count+=length
                        new_seq = ("N"*length) + new_seq[stop:]
                    else:
                        length = stop-start
                        new_seq = new_seq[:start] + ("N"*length) + new_seq[stop:]
                        n_count+=length

                n_diff = new_seq.count("N") - str(record.seq).count("N") 
                fw.write(f">{record.description}\n{new_seq}\n")
            else:
                fw.write(f">{record.description}\n{record.seq}\n")