#!/usr/bin/env python3

import os
import sys
import yaml
import collections
import csv

from piranha.utils import misc

from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *

def parse_barcodes_csv(barcodes_csv,config):

    misc.add_file_to_config(KEY_BARCODES_CSV,barcodes_csv,config)
    
    if not config[KEY_BARCODES_CSV]:
        sys.stderr.write(cyan(f"Error: No barcode csv file provided.\n"))
        sys.exit(-1)

    barcodes= []
    samples = []
    with open(config[KEY_BARCODES_CSV],"r") as f:
        reader = csv.DictReader(f)
        for col in ["barcode","sample"]:
            if not col in reader.fieldnames:
                sys.stderr.write(cyan(f"`{col}` must be a column name in barcode csv file.\n"))
                sys.exit(-1)

        for row in reader:
            if row["barcode"] in barcodes:
                barcode = row["barcode"]
                sys.stderr.write(cyan(f"`{barcode}` duplicated in barcode csv file. Note: barcodes must be unique.\n"))
                sys.exit(-1)
            if row["sample"] in samples:
                print(cyan(f"Warning: `{sample}` sample name provided for multiple barcodes."))
            barcodes.append(row["barcode"])

    config[KEY_BARCODES] = barcodes


def parse_read_dir(readdir,config):

    if readdir:
        readdir = os.path.abspath(readdir)

    misc.add_arg_to_config(KEY_READDIR,readdir,config)
    misc.check_path_exists(config[KEY_READDIR])

    run_id = False

    count_read_files = collections.Counter()
    for r,d,f in os.walk(config[KEY_READDIR]):
        for fn in f:
            if fn.endswith(".fastq"):
                if not run_id:
                    try:
                        # read run id from first file name it comes to
                        run_id = fn.split(".")[0].split("_")[-2]
                    except:
                        run_id = ""
                barcode = r.split("/")[-1]
                count_read_files[barcode]+=1
    
    config[KEY_RUNID] = run_id
    
    print(green("Found read files"))
    print(yellow("-----------------------"))
    for d in sorted(count_read_files):
        if count_read_files[d] == 1:
            print(green(f"Barcode {d}:\t") + f"{count_read_files[d]} fastq file")
        else:
            print(green(f"Barcode {d}:\t") + f"{count_read_files[d]} fastq files")

        if d not in config[KEY_BARCODES]:
            print(cyan(f"Warning: barcode {d} is not included in the input barcodes csv and will not be analysed."))

    for barcode in config[KEY_BARCODES]:
        if barcode not in count_read_files:
            sys.stderr.write(cyan(f"Error: No read files identified for barcode `{barcode}`.\nPlease check barcode csv has correct barcode names."))
            sys.exit(-1)

def parse_input_group(barcodes_csv,readdir,reference_sequences,config):

    parse_barcodes_csv(barcodes_csv,config)

    

    parse_read_dir(readdir,config)

    misc.add_file_to_config(KEY_REFERENCE_SEQUENCES,reference_sequences,config)
    misc.check_path_exists(config[KEY_REFERENCE_SEQUENCES])


