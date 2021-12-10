#!/usr/bin/env python3
from piranha.utils.log_colours import green,cyan
from piranha.utils import misc
from piranha.utils.config import *

import sys
import os
import glob

from datetime import datetime 
from datetime import date
import tempfile
import sys

"""
Desired behaviour

Default outdir -> analysis-2021-XX-YY,
check if that exists, append a number if it already exists -> analysis-2021-XX-YY_2
--overwrite flag to overwrite the previous output directory

If -o/--outdir then output that as the outdir instead, without a datestamp
If -p/--output-prefix then output directory as: prefix-2021-XX-YY

prefix-2021-XX-YY (or analysis-2021-XX-YY) will also be the name of the final analysis report.

"""

def datestamped_outdir(config):
    today = date.today()
    d = today.strftime("%Y-%m-%d")

    if not KEY_OUTDIR in config:
        expanded_path = os.path.expanduser(config[KEY_CWD])
        outdir = os.path.join(expanded_path, f"{config[KEY_OUTPUT_PREFIX]}_{d}")
    else:
        if config[KEY_DATESTAMP]:
            outdir = f"{config[KEY_OUTDIR]}_{d}"
        else:
            outdir = config[KEY_OUTDIR]

    if not config[KEY_OVERWRITE]:
        counter = 1
        while os.path.exists(outdir):
            if outdir.split("_")[-1].isdigit():
                outdir = "_".join(outdir.split("_")[:-1])
                outdir = f"{outdir}_{counter}"
            else:
                outdir = f"{outdir}_{counter}"
            counter +=1

    return outdir,d

def clear_old_files(config):
    if config[KEY_OVERWRITE] and os.path.exists(config[KEY_OUTDIR]):
        print(green("Overwriting previous output in ") + config[KEY_OUTDIR] + ".")
        old_files = glob.glob(f'{config[KEY_OUTDIR]}/*.*', recursive=True)

        for d in [old_files]:
            for f in d:
                try:
                    os.remove(f)
                    print(green("Removed: "),f)
                except:
                    print(cyan("Can't remove "),f)

def output_report_filename(d,config):
    if config[KEY_DATESTAMP]:
        output_report = f"{config[KEY_OUTPUT_PREFIX]}_{d}.html"
    else:
        output_report = f"{config[KEY_OUTPUT_PREFIX]}.html"
    return output_report

def set_up_tempdir(config):
    if config[KEY_NO_TEMP]:
        tempdir = config[KEY_OUTDIR]
        config[KEY_TEMPDIR] = tempdir
    elif KEY_TEMPDIR in config:
        to_be_dir = config[KEY_TEMPDIR]
        try:
            if not os.path.exists(to_be_dir):
                os.mkdir(to_be_dir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {to_be_dir}.\n'))
            sys.exit(-1)
        tempdir = tempfile.mkdtemp(dir=to_be_dir)
        config[KEY_TEMPDIR] = tempdir
    else:
        tempdir = tempfile.mkdtemp()
        config[KEY_TEMPDIR] = tempdir
        try:
            if not os.path.exists(tempdir):
                os.mkdir(tempdir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {tempdir}.\n'))
            sys.exit(-1)
        
        try:
            with open(os.path.join(tempdir, "test.txt"),"w") as fw:
                fw.write("Test")
        except:
            sys.stderr.write(cyan(f'Error: cannot write to temp directory {tempdir}.\n'))
            sys.exit(-1)

def output_group_parsing(outdir,output_prefix,overwrite,datestamp,tempdir,no_temp,config):
    
    misc.add_path_to_config(KEY_OUTDIR,outdir,config)
    misc.add_arg_to_config(KEY_OUTPUT_PREFIX,output_prefix,config)
    misc.add_arg_to_config(KEY_OVERWRITE,overwrite,config)
    misc.add_arg_to_config(KEY_DATESTAMP,datestamp,config)
    misc.add_path_to_config(KEY_TEMPDIR,tempdir,config)
    misc.add_arg_to_config(KEY_NO_TEMP,no_temp,config)

    config[KEY_OUTDIR],d = datestamped_outdir(config)
    
    set_up_tempdir(config)

    config[KEY_OUTPUT_REPORT] = output_report_filename(d,config)

    clear_old_files(config)

    if not os.path.exists(config[KEY_OUTDIR]):
        os.mkdir(config[KEY_OUTDIR])
