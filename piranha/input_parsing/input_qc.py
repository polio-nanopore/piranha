import os
import sys
import yaml
import collections
import csv

from pirahna.utils import misc

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

def parse_barcodes_csv(barcodes_csv,config):

    misc.add_file_to_config(KEY_BARCODES_CSV,barcodes_csv,config)
    
    barcodes = {}
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

            barcodes[row["barcode"]] = row["sample"]

    config[KEY_BARCODES] = barcodes


def parse_read_dir(readdir,config):
    misc.add_arg_to_config(KEY_READDIR,readdir,config)

    if readdir:
        readdir = os.path.abspath(readdir)

    misc.add_arg_to_config(KEY_READDIR,readdir,config)
    misc.check_path_exists(config[KEY_READDIR])

    count_read_files = {}
    for r,d,f in os.walk(config[KEY_READDIR]):
        for fn in f:
            if fn.endswith(".fastq"):
                count_read_files[d]+=1
    
    print(green("Found read files:"))
    for d in count_read_files:
        print(green(f"Barcode {d}:\t") + f"{count_read_files[d]} fastq files")
        if d not in config[KEY_BARCODES]:
            print(cyan(f"Warning: barcode {d} is not included in the input barcodes csv and will not be analysed."))

    for barcode in config[KEY_BARCODES]:
        if barcode not in count_read_files:
            sys.stderr.write(cyan(f"Error: No read files identified for barcode `{barcode}`.\nPlease check barcode csv has correct barcode names."))
            sys.exit(-1)

def parse_input_group(barcodes_csv,readdir,config):

    parse_barcodes_csv(barcodes_csv,config)
    parse_read_dir(readdir,config)


