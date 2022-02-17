#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv

from piranha.utils.config import *

def gather_fasta_files(summary_info, barcodes_csv, input_cns_list, output_file):
    analysis_info = collections.defaultdict(list)
    with open(summary_info, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            analysis_info[row[KEY_BARCODE]].append(row)

    input_metadata = {}
    with open(barcodes_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            input_metadata[row[KEY_BARCODE]] = row

    with open(output_file,"w") as fw:
        for cns_file in input_cns_list:
            for record in SeqIO.parse(cns_file, "fasta"):
                cns_info= record.description.split(" ")
                ref,barcode,var_count,var_string=cns_info[0].split("|")
                
                info = []
                for row in analysis_info[barcode]:
                    if row[KEY_REFERENCE] == ref:
                        info = row
                metadata = input_metadata[barcode]

                record_id = f"{metadata[KEY_SAMPLE]}|{barcode}|{info[KEY_REFERENCE_GROUP]}|{var_count}|{var_string}|{metadata[KEY_DATE]}"
                fw.write(f">{record_id}\n{record.seq}\n")

