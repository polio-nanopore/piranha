#!/usr/bin/env python3

import os
import sys
import yaml
import collections
import csv
from Bio import SeqIO

from piranha.utils import misc

from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *

def define_valid_wells():
    valid_wells = set()
    for i in range(1,13):
        for j in ["A","B","C","D","E","F","G","H"]:
            if i<10:
                well = f"{j}0{i}"
            else:
                well = f"{j}{i}"
            valid_wells.add(well)
    return valid_wells


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
            if row[KEY_BARCODE] and not row[KEY_SAMPLE]:
                continue
            if row[KEY_BARCODE] in barcodes:
                barcode = row[KEY_BARCODE]
                sys.stderr.write(cyan(f"`{barcode}` duplicated in barcode csv file. Note: barcodes must be unique.\n"))
                sys.exit(-1)
            if row[KEY_SAMPLE] in samples:
                print(cyan(f"Warning: `{row[KEY_SAMPLE]}` sample name provided for multiple barcodes."))
            for special_character in ["|",","," ",";"]:
                if special_character in row[KEY_BARCODE]:
                    sys.stderr.write(cyan(f"Special character `{special_character}` cannot be used in barcode or sample name. Please remove this character from barcode `{row[KEY_BARCODE]}` and restart.\n"))
                    sys.exit(-1)
                elif special_character in row[KEY_SAMPLE]:
                    sys.stderr.write(cyan(f"Special character `{special_character}` cannot be used in barcode or sample name. Please remove this character from sample `{row[KEY_SAMPLE]}` and restart.\n"))
                    sys.exit(-1)

            barcodes.append(row[KEY_BARCODE])
            samples.append(row[KEY_SAMPLE])

    config[KEY_BARCODES] = barcodes
    config[KEY_SAMPLES] = samples
    
    valid_wells = define_valid_wells()
    invalid_wells = set()
    with open(config[KEY_BARCODES_CSV],"r") as f:
        reader = csv.DictReader(f)
        if "well" in reader.fieldnames:
            well_counter = collections.Counter()
            for row in reader:
                well_counter[row["well"]]+=1
                if row["well"] not in valid_wells:
                    invalid_wells.add(row["well"])

            duplicates = set()
            for i in well_counter:
                if well_counter[i] > 1:
                    duplicates.add(i)

            if duplicates:
                sys.stderr.write(cyan(f"Duplicate well coordinates specified in barcodes.csv file.\n"))
                for i in duplicates:
                    print(f"- {i}")
                sys.exit(-1)

            if invalid_wells:
                sys.stderr.write(cyan(f"Invalid well coordinates specified in barcodes.csv file.\n"))
                for i in invalid_wells:
                    print(f"- {i}")
                sys.exit(-1)


def phylo_group_parsing(run_phylo_arg, supplementary_sequences_arg):

    misc.add_arg_to_config(KEY_RUN_PHYLO,run_phylo_arg,config)

    if config[KEY_RUN_PHYLO] not in [True, False]:
        sys.stderr.write(cyan(f"run_phylo argument must be either True/False if specified through the config file.\n"))
        sys.exit(-1)
    
    if config[KEY_RUN_PHYLO]:

        if supplementary_sequences_arg:
            misc.add_file_to_config(KEY_SUPPLEMENTARY_SEQUENCES,supplementary_sequences_arg,config)
        
        if config[KEY_SUPPLEMENTARY_SEQUENCES]:
            misc.check_path_exists(config[KEY_SUPPLEMENTARY_SEQUENCES])

        try:
            records = SeqIO.index(config[KEY_SUPPLEMENTARY_SEQUENCES],"fasta")
            print(green("Supplementary sequences:"), len(records), "sequences parsed.")
        except:
            sys.stderr.write(cyan(f"Failed to parse supplementary sequence file, check it is in FASTA format.\n"))
            sys.exit(-1)

def parse_read_dir(readdir,config):

    if readdir:
        readdir = os.path.abspath(readdir)

    misc.add_arg_to_config(KEY_READDIR,readdir,config)
    misc.check_path_exists(config[KEY_READDIR])

    run_id = False

    count_read_files = collections.Counter()
    for r,d,f in os.walk(config[KEY_READDIR]):
        for fn in f:
            if fn.endswith(".fastq") or fn.endswith(".fq") or fn.endswith(".gz") or fn.endswith(".gzip"):
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
            print(green(f"Barcode {barcode}:\t") + f"0 fastq files")
            print(cyan(f"Warning: No read files identified for barcode `{barcode}`.\nThis may be a negative control or a failed sample, but be aware it will not be analysed."))

def parse_input_group(barcodes_csv,readdir,reference_sequences,config):

    parse_barcodes_csv(barcodes_csv,config)

    parse_read_dir(readdir,config)

    misc.add_file_to_config(KEY_REFERENCE_SEQUENCES,reference_sequences,config)
    misc.check_path_exists(config[KEY_REFERENCE_SEQUENCES])

    #check they have unique identifiers
    seq_ids = collections.Counter()
    for record in SeqIO.parse(config[KEY_REFERENCE_SEQUENCES],"fasta"):
        seq_ids[record.id]+=1
    
    more_than_once = []
    for seq in seq_ids:
        if seq_ids[seq]>1:
            more_than_once.append(seq)
    print_str = "\n - ".join(more_than_once)
    if len(more_than_once)>0:
        sys.stderr.write(cyan(f"\nReference fasta file contains duplicate sequence IDs:\n"))
        sys.stderr.write(f" - ")
        sys.stderr.write(f"{print_str}\n")
        sys.stderr.write(cyan(f"Please remove duplicates from file and run again.\n"))
        sys.exit(-1)


def control_group_parsing(positive_control, negative_control, config):
    # mod to allow multiple pos and negative samples
    if positive_control:
        positive_control = positive_control.split(",")
    if negative_control:
        negative_control = negative_control.split(",")
        
    misc.add_arg_to_config(KEY_POSITIVE,positive_control,config)
    misc.add_arg_to_config(KEY_NEGATIVE,negative_control,config)

    for pos in config[KEY_POSITIVE]:
        if pos not in config[KEY_SAMPLES]:
            print(cyan(f"Warning: cannot find positive control in barcode csv file: {pos}"))
            
    for neg in config[KEY_NEGATIVE]:
        if neg not in config[KEY_SAMPLES]:
            print(cyan(f"Warning: cannot find negative control in  barcode csv file: {neg}"))





