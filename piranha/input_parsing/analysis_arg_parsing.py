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

def analysis_group_parsing(min_read_length,max_read_length,min_read_depth,min_read_pcent,config):

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_MIN_READ_LENGTH,min_read_length,config)
    misc.add_arg_to_config(KEY_MAX_READ_LENGTH,max_read_length,config)
    misc.add_arg_to_config(KEY_MIN_READS,min_read_depth,config)
    misc.add_arg_to_config(KEY_MIN_PCENT,min_read_pcent,config)

    for key in [KEY_MIN_READ_LENGTH,KEY_MAX_READ_LENGTH,KEY_MIN_READS,KEY_MIN_PCENT]:
        check_if_int(key,config)

def sample_type(sample_type_arg,config):
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_SAMPLE_TYPE,sample_type_arg,config)
    
    if config[KEY_SAMPLE_TYPE] not in valid_sample_types:
        error_str = ', '.join(valid_sample_types)
        sys.stderr.write(cyan(f"`{KEY_SAMPLE_TYPE}` must be one of {error_str}.\n"))
        sys.exit(-1)

def analysis_mode(analysis_mode_arg,config):
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_ANALYSIS_MODE,analysis_mode_arg,config)

    if config[KEY_ANALYSIS_MODE] not in valid_analysis_modes:
        error_str = ', '.join(valid_analysis_modes)
        sys.stderr.write(cyan(f"`{KEY_ANALYSIS_MODE}` must be one of {error_str}.\n"))
        sys.exit(-1)

    # known issue with this: if these values get overwritten by an input config file, the defaults will override them
    analysis_mode = config[KEY_ANALYSIS_MODE]
    if analysis_mode == VALUE_ANALYSIS_MODE_WG_2TILE and (config[KEY_MIN_READ_LENGTH]==READ_LENGTH_DEFAULT_VP1[0] and config[KEY_MIN_READ_LENGTH]==READ_LENGTH_DEFAULT_VP1[1]):
        config[KEY_MIN_READ_LENGTH] = READ_LENGTH_DEFAULT_WG_2TILE[0]
        config[KEY_MAX_READ_LENGTH] = READ_LENGTH_DEFAULT_WG_2TILE[1]

    # print(config[KEY_MIN_READ_LENGTH],config[KEY_MAX_READ_LENGTH])
    # if config[KEY_ANALYSIS_MODE] != "vp1":
    #     sys.stderr.write(cyan(f"Only `vp1` analysis mode currently implemented.\n"))
    #     sys.exit(-1)