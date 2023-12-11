#!/usr/bin/env python3

import os
import sys
import yaml
import subprocess
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

from piranha.utils import misc
import piranha.utils.custom_logger as custom_logger
import piranha.utils.log_handler_handle as lh

def check_if_int(key,config):
    if config[key]:
        try:
            config[key] = int(config[key])
        except:
            sys.stderr.write(cyan(f"`{key}` must be numerical.\n"))
            sys.exit(-1)

def check_if_float(key,config):
    if config[key]:
        try:
            config[key] = float(config[key])
        except:
            sys.stderr.write(cyan(f"`{key}` must be numerical.\n"))
            sys.exit(-1)

def get_available_medaka_models():
    models = []
    result = subprocess.run(["medaka", "tools", "list_models"],stdout=subprocess.PIPE)
    result_string = result.stdout.decode('utf-8')
    for line in result_string.split("\n"):
        if line.startswith("Available"):
            models = line.split(": ")[1].split(", ")
    return models

def medaka_options_parsing(medaka_model,medaka_list_models,config):
    models = get_available_medaka_models()

    if medaka_list_models:
        print(green("Available medaka models:"))
        for i in models:
            print(f"- {i}")
        sys.exit(0)

    misc.add_arg_to_config(KEY_MEDAKA_MODEL,medaka_model,config)
    if config[KEY_MEDAKA_MODEL] not in models:
        sys.stderr.write(cyan(f"Medaka model specified not valid: `{config[KEY_MEDAKA_MODEL]}`.\nPlease use --medaka-list-models to see which models are available.\nIf needed, update medaka version with `pip install --upgrade medaka`.\n"))
        sys.exit(-1)

def analysis_group_parsing(min_read_length,max_read_length,min_read_depth,min_read_pcent,min_aln_block,primer_length,min_map_quality,config):

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_MIN_READ_LENGTH,min_read_length,config)
    misc.add_arg_to_config(KEY_MAX_READ_LENGTH,max_read_length,config)
    misc.add_arg_to_config(KEY_MIN_READS,min_read_depth,config)
    misc.add_arg_to_config(KEY_MIN_PCENT,min_read_pcent,config)
    misc.add_arg_to_config(KEY_MIN_ALN_BLOCK,min_aln_block,config)
    misc.add_arg_to_config(KEY_PRIMER_LENGTH,primer_length,config)
    misc.add_arg_to_config(KEY_MIN_MAP_QUALITY,min_map_quality,config)

    for key in [KEY_MIN_READ_LENGTH,KEY_MAX_READ_LENGTH,KEY_MIN_READS,KEY_MIN_ALN_BLOCK]:
        check_if_int(key,config)
    
    check_if_float(KEY_MIN_PCENT,config)

def sample_type(sample_type_arg,config):
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_SAMPLE_TYPE,sample_type_arg,config)
    
    if config[KEY_SAMPLE_TYPE] not in VALID_SAMPLE_TYPES:
        error_str = ', '.join(VALID_SAMPLE_TYPES)
        sys.stderr.write(cyan(f"`{KEY_SAMPLE_TYPE}` must be one of {error_str}.\n"))
        sys.exit(-1)

def analysis_mode(analysis_mode_arg,config):
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_ANALYSIS_MODE,analysis_mode_arg,config)

    if config[KEY_ANALYSIS_MODE] not in VALID_ANALYSIS_MODES:
        error_str = ', '.join(VALID_ANALYSIS_MODES)
        sys.stderr.write(cyan(f"`{KEY_ANALYSIS_MODE}` must be one of {error_str}.\n"))
        sys.exit(-1)


    #must do this before parsing args for min and max read len
    if config[KEY_ANALYSIS_MODE] in READ_LENGTH_DICT:
        config[KEY_MIN_READ_LENGTH] = READ_LENGTH_DICT[config[KEY_ANALYSIS_MODE]][0]
        config[KEY_MAX_READ_LENGTH] = READ_LENGTH_DICT[config[KEY_ANALYSIS_MODE]][1]

    print(green(f"Default read length filter for {config[KEY_ANALYSIS_MODE]}:") + f" {config[KEY_MIN_READ_LENGTH]}-{config[KEY_MAX_READ_LENGTH]}")
    # if config[KEY_ANALYSIS_MODE] != "vp1":
    #     sys.stderr.write(cyan(f"Only `vp1` analysis mode currently implemented.\n"))
    #     sys.exit(-1)



def haplo_group_parsing(run_haplotyping,
                        haplotype_sample_size,
                        min_allele_frequency,
                        max_haplotypes,
                        min_haplotype_distance,
                        min_haplotype_depth,
                        config):

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_RUN_HAPLOTYPING,run_haplotyping,config)
    misc.add_arg_to_config(KEY_HAPLOTYPE_SAMPLE_SIZE,haplotype_sample_size,config)
    misc.add_arg_to_config(KEY_MIN_ALLELE_FREQUENCY,min_allele_frequency,config)
    misc.add_arg_to_config(KEY_MAX_HAPLOTYPES,max_haplotypes,config)
    misc.add_arg_to_config(KEY_MIN_HAPLOTYPE_DISTANCE,min_haplotype_distance,config)
    misc.add_arg_to_config(KEY_MIN_HAPLOTYPE_DEPTH,min_haplotype_depth,config)

    for key in [KEY_HAPLOTYPE_SAMPLE_SIZE,KEY_MAX_HAPLOTYPES,KEY_MIN_HAPLOTYPE_DISTANCE,KEY_MIN_HAPLOTYPE_DEPTH]:
        check_if_int(key,config)
    
    check_if_float(KEY_MIN_ALLELE_FREQUENCY,config)
