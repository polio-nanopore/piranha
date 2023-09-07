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
    
