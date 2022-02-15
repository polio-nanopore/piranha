import os
import sys
import yaml
from datetime import date

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

from piranha.utils import misc
import piranha.utils.custom_logger as custom_logger
import piranha.utils.log_handler_handle as lh


def get_defaults():
    today = date.today()
    default_dict = {
                    
                    #input options
                    KEY_READDIR: False,
                    KEY_REFERENCE_SEQUENCES: False,
                    KEY_BARCODES_CSV: False,
                    KEY_RUNID: False,

                    # Output options
                    KEY_OUTPUT_PREFIX:VALUE_OUTPUT_PREFIX,
                    KEY_DATESTAMP:False,
                    KEY_NO_TEMP:False,
                    KEY_OVERWRITE:False,
                    KEY_REFERENCES_FOR_CNS:VALUE_REFERENCES_FOR_CNS,
                    KEY_SUMMARY_HEADERS: VALUE_SUMMARY_HEADERS,

                    # input seq options 
                    KEY_ANALYSIS_MODE:VALUE_ANALYSIS_MODE, #options are stool or environmental
                    KEY_MIN_READ_LENGTH:VALUE_MIN_READ_LENGTH,
                    KEY_MAX_READ_LENGTH:VALUE_MAX_READ_LENGTH,
                    KEY_MIN_READS:VALUE_MIN_READS,   # where to pad to using datafunk
                    KEY_MIN_PCENT:VALUE_MIN_PCENT,

                    # misc defaults
                    KEY_THREADS:1,
                    KEY_VERBOSE:False,

                    KEY_COLOUR_MAP: VALUE_COLOUR_MAP,
                    KEY_COLOUR_THEME: VALUE_COLOUR_THEME

                    }
    return default_dict

def valid_args():
    return [
        KEY_READDIR,
        KEY_BARCODES_CSV,
        KEY_OUTDIR,
        KEY_OUTPUT_PREFIX,
        KEY_TEMPDIR,
        KEY_NO_TEMP,
        KEY_OVERWRITE,
        KEY_MIN_READ_LENGTH,
        KEY_MAX_READ_LENGTH,
        KEY_MIN_READS,
        KEY_MIN_PCENT,
        KEY_THREADS,
        KEY_VERBOSE,
        KEY_REFERENCE_SEQUENCES,
        KEY_REFERENCES_FOR_CNS
    ]

def check_configfile(cwd,config_arg):
    configfile = os.path.join(cwd,config_arg)

    ending = configfile.split(".")[-1]

    if ending not in ["yaml","yml"]:
        sys.stderr.write(cyan(f'Error: config file {configfile} must be in yaml format.\n'))
        sys.exit(-1)
    
    elif not os.path.isfile(configfile):
        sys.stderr.write(cyan(f'Error: cannot find config file at {configfile}\n'))
        sys.exit(-1)
    else:
        print(green(f"Input config file:") + f" {configfile}")
        return configfile

def load_yaml(f):
    try:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
    except:
        sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
        sys.exit(-1)
    return input_config

def return_path_keys():
    return [KEY_READDIR,KEY_OUTDIR,KEY_TEMPDIR]

def setup_absolute_paths(path_to_file,value):
    return os.path.join(path_to_file,value)


def parse_yaml_file(configfile,configdict):
    overwriting = 0
    path_keys = return_path_keys()

    path_to_file = os.path.abspath(os.path.dirname(configfile))
    configdict[KEY_INPUT_PATH] = path_to_file
    valid_keys = valid_args()

    invalid_keys = []
    with open(configfile,"r") as f:
        input_config = load_yaml(f) # try load file else exit with msg

        for key in input_config:
            value = input_config[key]
            if value == None: # dont count blank entries
                pass
            else:
                clean_key = key.lstrip("-").replace("-","_").rstrip(" ").lstrip(" ").lower()

                if clean_key not in valid_keys:
                    invalid_keys.append(key)
                    break
                    
                if clean_key in path_keys:
                    value = setup_absolute_paths(path_to_file,value)

                configdict[clean_key] = value
                overwriting += 1

    if len(invalid_keys)==1:
        sys.stderr.write(cyan(f'Error: invalid key in config file.\n') + f'\t- {invalid_keys[0]}\n')
        sys.exit(-1)
    elif len(invalid_keys) >1:
        keys = ""
        for i in invalid_keys:
            keys += f"\t- {i}\n"
        sys.stderr.write(cyan(f'Error: invalid keys in config file.\n') + f'{keys}')
        sys.exit(-1)
    print(green(f"Adding {overwriting} arguments to internal config."))

def setup_config_dict(cwd,config_arg):
    config = get_defaults()
    config[KEY_CWD] = cwd
    if config_arg:
        configfile = check_configfile(cwd,config_arg)
        parse_yaml_file(configfile,config)

    else:
        config[KEY_INPUT_PATH] = cwd
    return config

def misc_args_to_config(verbose,threads, config):
    misc.add_arg_to_config(KEY_VERBOSE,verbose,config)
    misc.add_arg_to_config(KEY_THREADS,threads,config)

def set_up_verbosity(config):
    if config[KEY_VERBOSE]:
        config[KEY_QUIET] = False
        config[KEY_LOG_API] = ""
        config[KEY_LOG_STRING] = ""
    else:
        config[KEY_QUIET] = True
        config[KEY_LOG_API] = ""
        lh_path = os.path.realpath(lh.__file__)
        config[KEY_LOG_STRING] = f"--quiet --log-handler-script {lh_path} "