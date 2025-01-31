#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt

import snakemake


import piranha.utils.custom_logger as custom_logger
from piranha.utils.log_colours import green,cyan,red
from piranha.utils.config import *

def run_snakemake(snake_config,snakefile,config):
    if config[KEY_VERBOSE]:
        print(red("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR], config=snake_config, cores=config[KEY_THREADS],lock=False
                                    )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR], config=snake_config, cores=config[KEY_THREADS],lock=False,
                                    quiet=True,log_handler=logger.log_handler
                                    )
    return status

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

def add_check_valid_arg(KEY,arg,valid_values,config):
    add_arg_to_config(KEY,arg,config)
    if config[KEY] not in valid_values:
        sys.stderr.write(cyan(f'{config[KEY]} not a valid option for {KEY}. Please specify one of `{"`, `".join(valid_values)}`\n'))
        sys.exit(-1)

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
    elif arg == 0:
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

def check_path_exists(path):
    if not os.path.exists(path):
        sys.stderr.write(cyan(f"Error: path {path} does not exist. Check input.\n"))
        sys.exit(-1)



def header(v):
    print(green("""\n
                            __                      __            
                    ______ |__|_____ _____  |\____ |  |__ _____   
                    \____ \|  |\  _ \\\__  \ |     \|  |  \\\__  \  
                    |  |_| |  ||  |\/ / __ \|   |  \   |  \/ __ \ 
                    |   __/|__||__|  (____ / ___|__/___|___(____/ 
                    |__|                                          

                     **** Poliovirus Investigation Resource ****
                   **** Automating Nanopore Haplotype Analysis ****
                """)+green(f"""
                                        {v}""")+green("""
                        ****************************************""")+("""
                             Developed by researchers at""")+green(f"""
                               University of Edinburgh""")+("""
                                in collaboration with""")+green(f"""
                           Imperial College London & NIBSC""")+("""
                                      Supported by:""")+green(f"""
                    ARTIC Network Wellcome Trust Collaborators Award
                                    (206298/Z/17/Z)
                           Bill and Melinda Gates Foundation 
                                    (OPP1207299)

\n"""))

def preamble(v):
    header(v)

def minimap2_help():
    print(green("""\n
*** minimap2 [options] within piranha ***
Options:
  Indexing:
    -k INT       k-mer size (no larger than 28) [15]
    -w INT       minimizer window size [10]
  Mapping:
    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [0.0002]
    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]
    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]
    -r NUM       bandwidth used in chaining and DP-based alignment [500]
    -n INT       minimal number of minimizers on a chain [3]
    -m INT       minimal chaining score (matching bases minus log gap penalty) [40]
  Alignment:
    -A INT       matching score [2]
    -B INT       mismatch penalty [4]
    -O INT[,INT] gap open penalty [4,24]
    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
    -z INT[,INT] Z-drop score and inversion Z-drop score [400,200]
    -s INT       minimal peak DP alignment score [80]
    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
  Preset:
    -x STR       preset (always applied before other options; see minimap2.1 for details) []
                 - map-pb/map-ont: PacBio/Nanopore vs reference mapping
                 - ava-pb/ava-ont: PacBio/Nanopore read overlap
                 - asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
                 - splice: long-read spliced alignment
                 - sr: genomic short-read mapping

See `man ./minimap2.1' for detailed description of these and other advanced command-line options.

\n"""))

