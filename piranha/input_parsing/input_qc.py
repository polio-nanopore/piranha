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

# def qc_supplementary_sequence_file(supplementary_sequences):
#     misc.check_path_exists(supplementary_sequences)
#     incorrect = 0
#     total = 0
#     seq_ids = set()
#     try:
#         for record in SeqIO.parse(supplementary_sequences,"fasta"):
#             seq_ids.add(record.id)

#             total +=1
#             passed = False
#             for field in record.description.split(" "):
#                 if field.startswith(KEY_DISPLAY_NAME):
#                     passed=True
#                     continue
#             if not passed:
#                 incorrect +=1
#     except:
#         sys.stderr.write(cyan(f"Failed to parse supplementary sequence file, check it is in FASTA format.\n"))
#         sys.exit(-1)

#     if incorrect >= 1:
#         sys.stderr.write(cyan(f"Supplementary sequences file lacks `{KEY_DISPLAY_NAME}` annotation in header of {incorrect} out of {total} sequences parsed.\n"))
#         sys.exit(-1)
#     else:
#         print(green("Supplementary sequences:"), total, "sequences parsed.")
    
#     return seq_ids

def qc_supplementary_metadata_file(supplementary_metadata,seq_ids,config):

    with open(supplementary_metadata,"r") as f:
        reader = csv.DictReader(f)

        missing = set()
        for col in config[KEY_SUPPLEMENTARY_METADATA_COLUMNS]:
            if col not in reader.fieldnames and col not in VALUE_SUPPLEMENTARY_METADATA_COLUMNS:
                missing.add(col)
        if missing:
            sys.stderr.write(cyan(f"Error: {supplementary_metadata} missing specified columns:\n"))
            for i in missing:
                sys.stderr.write(cyan(f"- {i}\n"))
            sys.exit(-1)


def parse_fasta_file(supplementary_datadir,supp_file,seq_records,no_reference_group,total_seqs,seq_info,config):
    for record in SeqIO.parse(os.path.join(supplementary_datadir,supp_file),"fasta"):
        total_seqs["total"] +=1
        ref_group = ""
        for field in record.description.split(" "):
            if field.startswith(KEY_DISPLAY_NAME):
                ref_group = field.split("=")[1]
        
        if ref_group not in config[KEY_REFERENCES_FOR_CNS]:
            no_reference_group.add(record.id)
        else:
            total_seqs[reference_group]+=1
            seq_records.append(record)
            seq_info[record.id]= {}

def check_there_are_seqs(total_seqs,supplementary_datadir,no_reference_group,config):
    if total_seqs["total"]==0:
        sys.stderr.write(cyan(f"Error: No sequence files matched in `{supplementary_datadir}`.\nEnsure the directory provided contains FASTA files with appropriate annotations in the header.\n"))
        sys.stderr.write(cyan(f"Header must specify one of {config[KEY_REFERENCES_FOR_CNS]} under {KEY_DISPLAY_NAME}=X, where X is the appropriate group to be included in phylo pipeline.\n"))
        sys.exit(-1)

    elif no_reference_group:
        print(cyan(f"Warning: not all sequences in {supplementary_datadir} included in analysis because they lack annotations in the header description."))
        print(cyan("Number of seqs not included:"),f"{len(no_reference_group)}")

    print(green(f"Total supplementary sequences:"),total_seqs["total"])
    for reference_group in config[KEY_REFERENCES_FOR_CNS]:
        print(green(f"Total {reference_group} sequences:"), total_seqs[reference_group])

def gather_supplementary_data(supplementary_datadir,supplementary_sequences,supplementary_metadata,supplementary_metadata_id_column,config):

    seq_records = []
    
    no_reference_group = set()
    total_seqs = collections.Counter()

    seq_info = collections.defaultdict(dict)
    metadata_info = collections.defaultdict(dict)
    passed_over_csvs = set()
    with open(supplementary_sequences,"w") as fw:
        for r,d,f in os.walk(supplementary_datadir):
            for supp_file in f:
                if supp_file.endswith(".fasta") or supp_file.endswith(".fa"):
                    parse_fasta_file(supplementary_datadir,supp_file,seq_records,no_reference_group,total_seqs,seq_info,config)
                
                elif supp_file.endswith("csv"):
                    with open(os.path.join(supplementary_datadir,supp_file),"r")  as f:
                        reader = csv.DictReader(f)
                        if supplementary_metadata_id_column not in reader.fieldnames:
                            passed_over_csvs.add(supp_file)
                            continue
                        else:
                            for row in reader:
                                seq_name = row[supplementary_metadata_id_column]
                                metadata_info[seq_name] = row
        
        check_there_are_seqs(total_seqs,supplementary_datadir,no_reference_group,config)

        SeqIO.write(seq_records,fw, "fasta")

        supplementary_metadata_header = set()

        for seq_id in seq_info:
            if seq_id in metadata_info:
                row = metadata_info[seq_id]
                seq_info[seq_id] = row
                for col in row:
                    supplementary_metadata_header.add(col)

        with open(supplementary_metadata,"w") as f:
            writer = csv.DictWriter(fw, fieldnames=supplementary_metadata_header,lineterminator="\n")
            writer.writeheader()
            for seq_id in seq_info:
                row = seq_info[seq_id]
                for field in supplementary_metadata_header:
                    if field not in row:
                        row[field] = ""
                
                writer.writerow(row)

def phylo_group_parsing(run_phylo_arg,
                        update_local_database,
                        supplementary_datadir,
                        phylo_metadata_columns_arg,
                        barcodes_csv,
                        supplementary_metadata_columns_arg,
                        supplementary_metadata_id_column_arg,
                        config):

    misc.add_arg_to_config(KEY_RUN_PHYLO,run_phylo_arg,config)

    if config[KEY_RUN_PHYLO] not in [True, False]:
        sys.stderr.write(cyan(f"`{KEY_RUN_PHYLO}` argument must be either True/False if specified through the config file.\n"))
        sys.exit(-1)
    
    if config[KEY_RUN_PHYLO]:

        misc.add_path_to_config(KEY_SUPPLEMENTARY_DATADIR,supplementary_datadir,config)
        misc.check_path_exists(config[KEY_SUPPLEMENTARY_DATADIR])

        misc.add_arg_to_config(KEY_UPDATE_LOCAL_DATABASE,update_local_database,config)
        if config[KEY_UPDATE_LOCAL_DATABASE] not in [True, False]:
            sys.stderr.write(cyan(f"`{KEY_UPDATE_LOCAL_DATABASE}` argument must be either True/False if specified through the config file.\n"))
            sys.exit(-1)

        if config[KEY_UPDATE_LOCAL_DATABASE] and not config[KEY_SUPPLEMENTARY_DATADIR]:
            sys.stderr.write(cyan(f"Error: Cannot update local database with new sequences as no supplementary database has been provided.\n"))
            sys.exit(-1)

        misc.add_arg_to_config(KEY_SUPPLEMENTARY_METADATA_ID_COLUMN,supplementary_metadata_id_column_arg,config)
        
        if config[KEY_SUPPLEMENTARY_DATADIR]:
            db_dir = os.path.join(config[KEY_TEMPDIR],"local_db")
            if not os.path.exists(db_dir):
                os.mkdir(db_dir)

            config[KEY_SUPPLEMENTARY_SEQUENCES] = os.path.join(config[KEY_TEMPDIR],"local_db","supp_seqs.fasta")
            config[KEY_SUPPLEMENTARY_METADATA] = os.path.join(config[KEY_TEMPDIR],"local_db","supp_metadata.csv")
            
            gather_supplementary_data(config[KEY_SUPPLEMENTARY_DATADIR],config[KEY_SUPPLEMENTARY_SEQUENCES],config[KEY_SUPPLEMENTARY_METADATA],config[KEY_SUPPLEMENTARY_METADATA_ID_COLUMN],config)

        seq_ids = set()
        
        misc.add_arg_to_config(KEY_SUPPLEMENTARY_METADATA_COLUMNS,supplementary_metadata_columns_arg,config)
        if not type(config[KEY_SUPPLEMENTARY_METADATA_COLUMNS])==list:
            config[KEY_SUPPLEMENTARY_METADATA_COLUMNS] = [config[KEY_SUPPLEMENTARY_METADATA_COLUMNS]]

        if config[KEY_SUPPLEMENTARY_METADATA]:
            qc_supplementary_metadata_file(config[KEY_SUPPLEMENTARY_METADATA],
                                            seq_ids,
                                            config)

        misc.add_arg_to_config(KEY_PHYLO_METADATA_COLUMNS,phylo_metadata_columns_arg,config)
        with open(barcodes_csv,"r") as f:
            reader = csv.DictReader(f)
            missing = set()
            for col in config[KEY_PHYLO_METADATA_COLUMNS]:
                if col not in reader.fieldnames and col not in VALUE_PHYLO_METADATA_COLUMNS:
                    missing.add(col)

            if missing:
                sys.stderr.write(cyan(f"The following {KEY_PHYLO_METADATA_COLUMNS} columns are missing from the barcodes.csv file:\n"))
                for i in missing:
                    sys.stderr.write(cyan(f"- {i}\n"))
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

    misc.add_arg_to_config(KEY_POSITIVE,positive_control,config)
    misc.add_arg_to_config(KEY_NEGATIVE,negative_control,config)

    if config[KEY_POSITIVE] and not type(config[KEY_POSITIVE]) == list:
        config[KEY_POSITIVE] = config[KEY_POSITIVE].split(",")
    if config[KEY_NEGATIVE] and not type(config[KEY_NEGATIVE]) == list:
        config[KEY_NEGATIVE] = config[KEY_NEGATIVE].split(",")

    for pos in config[KEY_POSITIVE]:
        if pos not in config[KEY_SAMPLES]:
            print(cyan(f"Warning: cannot find positive control in barcode csv file: {pos}"))
            
    for neg in config[KEY_NEGATIVE]:
        if neg not in config[KEY_SAMPLES]:
            print(cyan(f"Warning: cannot find negative control in  barcode csv file: {neg}"))





