#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt

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



def add_col_to_metadata(new_column_name, new_column_dict, metadata, new_metadata, match_column, config): 
    #dictionary currently is key=sequence name and value=new col value
    print(green("Adding column to master metadata:"), new_column_name)
    with open(new_metadata, 'w') as fw:
        
        with open(metadata,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            header.append(new_column_name)
            config[KEY_QUERY_CSV_HEADER] = header
            
            writer = csv.DictWriter(fw, fieldnames=config[KEY_QUERY_CSV_HEADER],lineterminator='\n')
            writer.writeheader()

            for row in reader:
                new_row = row

                if new_row[match_column] in new_column_dict:
                    new_row[new_column_name] = new_column_dict[new_row[match_column]]
                else:
                    new_row[new_column_name] = ""

                writer.writerow(new_row)

def check_date_format(to_check, line_count, header):

    date_format = "%Y-%m-%d"
    
    try:
        dt.datetime.strptime(to_check, date_format).date()
    except:
        sys.stderr.write(cyan(f'Date {to_check} on line {line_count} in column {header} in incorrect format. Please use YYYY-MM-DD\n'))
        sys.exit(-1)
                
def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

def add_file_to_config(key,arg,config):
    if arg:
        path_to_file = os.path.abspath(config["cwd"])
        full_path = os.path.join(path_to_file,arg)
        config[key]=full_path

def add_path_to_config(key,arg,config):
    if arg:
        expanded_path = os.path.expanduser(arg)
        path_to_cwd = os.path.abspath(config["cwd"])
        full_path = os.path.join(path_to_cwd,expanded_path)
        config[key]=full_path
