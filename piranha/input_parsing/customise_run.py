#!/usr/bin/env python3

import os

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

def look_for_basecalled_reads(read_path_arg,cwd,config):

    if read_path_arg:
        expanded_path = os.path.expanduser(read_path_arg)
        read_path = os.path.join(cwd, expanded_path)
        if not os.path.exists(read_path):
            sys.stderr.write(cyan('Error: cannot find reads at {}\n'.format(read_path)))
            sys.exit(-1)
        else:
            config["read_path"] = read_path

    elif "read_path" in config:
        if config["read_path"]:
            expanded_path = os.path.expanduser(config["read_path"])
            read_path = os.path.join(config["path_to_config"], expanded_path)
            config["read_path"] = read_path
            fq_files = 0
            for r,d,f in os.walk(read_path):
                for fn in f:
                    filename = fn.lower()
                    if filename.endswith(".fastq") or filename.endswith(".fq"):
                        fq_files +=1

            if fq_files > 0:
                print(f"Found {fq_files} fastq files in the input directory")
            else:
                sys.stderr.write(cyan('Error: cannot find fastq files at {}\nPlease check your `--read-dir`'.format(read_path)))
                sys.exit(-1)
    else:
        sys.stderr.write(cyan('Error: `--read-dir` needed. Please input the path to the fastq read files either in the config file or via the command line.\n'))
        sys.exit(-1)

def look_for_barcodes_csv(barcodes_csv_arg,cwd,config):
    barcodes_csv = ""
    if barcodes_csv_arg:
    
        barcodes_csv = os.path.join(cwd,barcodes_csv_arg)
        if not os.path.exists(barcodes_csv):
            sys.stderr.write('Error: cannot find barcodes csv at {}\n'.format(barcodes_csv))
            sys.exit(-1)
    elif "barcodes_csv" in config:
        if config["barcodes_csv"]:
            expanded_path = os.path.expanduser(config["barcodes_csv"])
            barcodes_csv = os.path.join(config["path_to_config"], expanded_path)
    
    if barcodes_csv:
        print(f"Input barcodes csv file: {barcodes_csv}")
        barcodes = []
        with open(barcodes_csv, newline="") as f:
            reader = csv.DictReader(f)
            column_names = reader.fieldnames
            if "barcode" not in column_names:
                sys.stderr.write(f"Error: Barcode file missing header field `barcode`\n")
                sys.exit(-1)
            for row in reader: 
                if row["barcode"].startswith("NB") or row["barcode"].startswith("BC"):
                    barcodes.append(row["barcode"])
                else:
                    sys.stderr.write(f"Error: Please provide barcodes in the format `NB01` or `BC01`\n")
                    sys.exit(-1)
                    
            print(f"{len(barcodes)} barcodes read in from file")
            for i in barcodes:
                print(f"  - {i}")
            barcodes = ",".join(barcodes)
            config['barcodes_csv'] = barcodes_csv
            config["barcodes"] = barcodes
    else:
        config['barcodes_csv'] = ""
        config["barcodes"] = ""
        print(green(f"Note: No barcodes csv input"))

def check_barcode_kit():

    add_arg_to_config("barcode_kit",barcode_kit_arg, config)
    barcode_kit = config["barcode_kit"]
    if args.barcode_kit.lower() in ["native","pcr","rapid","all"]:
        config["barcode_set"] = args.barcode_kit.lower()
    else:
        sys.stderr.write(f"Error: Please enter a valid barcode kit: one of\n\t-native\n\t-pcr\n\t-rapid\n\t-all\n")
        sys.exit(-1)
