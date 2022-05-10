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

    if not os.path.exists(config[KEY_BARCODES_CSV]):
        sys.stderr.write(cyan(f"Error: Cannot find input file {config[KEY_BARCODES_CSV]}.\n"))
        sys.exit(-1)

    barcodes= []
    samples = []
    with open(config[KEY_BARCODES_CSV],"r") as f:
        reader = csv.DictReader(f)
        for col in [KEY_BARCODE,KEY_SAMPLE]:
            if not col in reader.fieldnames:
                sys.stderr.write(cyan(f"`{col}` must be a column name in barcode csv file.\n"))
                sys.exit(-1)

        for row in reader:
            if row[KEY_BARCODE] in barcodes:
                barcode = row[KEY_BARCODE]
                sys.stderr.write(cyan(f"`{barcode}` duplicated in barcode csv file. Note: barcodes must be unique.\n"))
                sys.exit(-1)
            if row[KEY_SAMPLE] in samples:
                print(cyan(f"Warning: `{sample}` sample name provided for multiple barcodes."))
            for special_character in ["|",","," ",";"]:
                if special_character in row[KEY_BARCODE]:
                    sys.stderr.write(cyan(f"`{special_character}` cannot be used in barcode or sample name. Please remove this character from barcode `{row[KEY_BARCODE]}` and restart.\n"))
                    sys.exit(-1)
                elif special_character in row[KEY_SAMPLE]:
                    sys.stderr.write(cyan(f"`{special_character}` cannot be used in barcode or sample name. Please remove this character from sample `{row[KEY_SAMPLE]}` and restart.\n"))
                    sys.exit(-1)
            barcodes.append(row[KEY_BARCODE])
            samples.append(row[KEY_SAMPLE])

    config[KEY_BARCODES] = barcodes
    config[KEY_SAMPLES] = samples


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

def control_group_parsing(positive_control, negative_control, config):
    misc.add_arg_to_config(KEY_POSITIVE,positive_control,config)
    misc.add_arg_to_config(KEY_NEGATIVE,negative_control,config)

    if config[KEY_POSITIVE] not in config[KEY_SAMPLES]:
        print(cyan(f"Warning: cannot find positive control: {config[KEY_POSITIVE]}"))
    if config[KEY_NEGATIVE] not in config[KEY_SAMPLES]:
        print(cyan(f"Warning: cannot find negative control: {config[KEY_NEGATIVE]}"))





