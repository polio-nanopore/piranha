#!/usr/bin/env python3
import subprocess
import os
import sys
from piranha.utils.log_colours import green,cyan
import importlib

def which(dependency):
    try:
        subprocess.check_output(["which", dependency])
        return True
    except subprocess.CalledProcessError:
        return False

def check_module(module, missing):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)

def check_this_dependency(dependency,missing):
    check = which(dependency)

    if not check:
        missing.append(dependency)


def look_for_guppy_barcoder(demultiplex_arg,path_to_guppy_arg,cwd,config):
    add_arg_to_config("demultiplex", demultiplex_arg, config)
    add_arg_to_config("path_to_guppy", path_to_guppy_arg, config)
    
    if config["demultiplex"]:
        if config["path_to_guppy"]:
            expanded_path = os.path.expanduser(config["path_to_guppy"])
            if config["path_to_guppy"].endswith("guppy_barcoder"):
                path_to_guppy = os.path.join(cwd,expanded_path)
            else:
                path_to_guppy = os.path.join(cwd,expanded_path,"guppy_barcoder")
                config["path_to_guppy"] = path_to_guppy
            os_cmd = os.system(f"{path_to_guppy} -v")

            if os_cmd != 0:
                sys.stderr.write(cyan(f'Error: guppy_barcoder at {path_to_guppy} fails to run\n'))
                sys.exit(-1)
        
        else:
            os_cmd = os.system(f"guppy_barcoder -v")
            if os_cmd != 0:
                sys.stderr.write(cyan(f'Error: please provide the path to guppy_barcoder (`--path-to-guppy`), add guppy_barcoder to your path, or run demultiplexing in MinKNOW\n'))
                sys.exit(-1)
            else:
                config["path_to_guppy"]


def check_dependencies(dependency_list, module_list):

    missing = []

    for dependency in dependency_list:
        check_this_dependency(dependency, missing)

    for module in module_list:
        check_module(module, missing)

    if missing:
        if len(missing)==1:
            sys.stderr.write(cyan(f'Error: Missing dependency `{missing[0]}`.')+'\nPlease update your piranha environment.\n')
            sys.exit(-1)
        else:
            dependencies = ""
            for i in missing:
                dependencies+=f"\t- {i}\n"

            sys.stderr.write(cyan(f'Error: Missing dependencies.')+f'\n{dependencies}Please update your piranha environment.\n')
            sys.exit(-1)
    else:
        print(green("All dependencies satisfied."))

# check_dependencies()
