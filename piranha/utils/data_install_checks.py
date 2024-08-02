#!/usr/bin/env python3
import pkg_resources
from piranha.utils.log_colours import green,cyan
import sys
import os
from piranha.utils.config import *
from piranha.utils import misc


def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('piranha', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(cyan(f'Error: Missing package data.')+f'\n\t- {filename}\nPlease install the latest piranha version with `piranha --update`.\n')
        sys.exit(-1)

def get_snakefile(thisdir,analysis_mode):

    snakefile = os.path.join(thisdir, 'scripts',f'piranha_{analysis_mode}.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n'))
        sys.exit(-1)
    return snakefile

def get_references(config):
    if config[KEY_ANALYSIS_MODE] == "vp1":
        package_data_check(REFERENCE_SEQUENCES_FILE_VP1,"data",KEY_REFERENCE_SEQUENCES,config)
    else:
        package_data_check(REFERENCE_SEQUENCES_FILE_WG,"data",KEY_REFERENCE_SEQUENCES,config)


def get_outgroups(config):

    if config[KEY_ANALYSIS_MODE] == "vp1":
        package_data_check(OUTGROUP_SEQUENCES_FILE_VP1,"data",KEY_OUTGROUP_SEQUENCES,config)
    else:
        package_data_check(OUTGROUP_SEQUENCES_FILE_WG,"data",KEY_OUTGROUP_SEQUENCES,config)

def check_install(language,config):
    misc.add_arg_to_config(KEY_LANGUAGE,language, config)
    if config[KEY_LANGUAGE] == "English":
        for resource in ENGLISH_RESOURCES:
            package_data_check(resource[RESOURCE_KEY_FILENAME],resource[RESOURCE_KEY_DIRECTORY],resource[RESOURCE_KEY],config)
    elif config[KEY_LANGUAGE] == "French":
        for resource in FRENCH_RESOURCES:
            package_data_check(resource[RESOURCE_KEY_FILENAME],resource[RESOURCE_KEY_DIRECTORY],resource[RESOURCE_KEY],config)
    else:
        sys.stderr.write(cyan(f'Error: `{config[KEY_LANGUAGE]}`not a valid language option. Available languages are English and French.\n'))
        sys.exit(-1)
    
    get_references(config)
    get_outgroups(config)


# config={}
# check_install()
