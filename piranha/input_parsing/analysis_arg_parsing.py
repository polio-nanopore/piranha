#!/usr/bin/env python3

import os
import sys
import yaml

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

def analysis_group_parsing(min_read_length,max_read_length,min_read_depth,min_read_pcent,config):

    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_MIN_READ_LENGTH,min_read_length,config)
    misc.add_arg_to_config(KEY_MAX_READ_LENGTH,max_read_length,config)
    misc.add_arg_to_config(KEY_MIN_READS,min_read_depth,config)
    misc.add_arg_to_config(KEY_MIN_PCENT,min_read_pcent,config)

    for key in [KEY_MIN_READ_LENGTH,KEY_MAX_READ_LENGTH,KEY_MIN_READS,KEY_MIN_PCENT]:
        check_if_int(key,config)

def analysis_mode(analysis_mode_arg,config):
    # if command line arg, overwrite config value
    misc.add_arg_to_config(KEY_ANALYSIS_MODE,analysis_mode_arg,config)

    if config[KEY_ANALYSIS_MODE] not in valid_analysis_modes:
        error_str = ', '.join(valid_analysis_modes)
        sys.stderr.write(cyan(f"`{KEY_ANALYSIS_MODE}` must be one of {error_str}.\n"))
        sys.exit(-1)

    if config[KEY_ANALYSIS_MODE] != "stool":
        sys.stderr.write(cyan(f"Only `stool` analysis mode currently implemented.\n"))
        sys.exit(-1)