import os
import sys
import yaml
import collections
import csv

from pirahna.utils import misc

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

def parse_barcodes_csv(barcodes_csv,config):


def parse_read_dir(readdir,config):
    misc.add_arg_to_config(KEY_READDIR,readdir,config)

    if readdir:
        readdir = os.path.abspath(readdir)

    misc.add_arg_to_config(KEY_READDIR,readdir,config)
    misc.check_path_exists(config[KEY_READDIR])

    count_read_files = {}
    for r,d,f in os.walk(config[KEY_READDIR]):
        for fn in f:
            if fn.endswith(".fastq"):
                count_read_files[d]+=1
    
    print(green("Found read files:"))
    for d in count_read_files:
        print(green(f"Barcode {d}:\t") + f"{count_read_files[d]} fastq files")


