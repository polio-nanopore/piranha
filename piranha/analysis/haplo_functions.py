#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

from piranha.utils.log_colours import green,cyan,red


def parse_partition_file(partition_file):

    partitions = collections.defaultdict(set)
    with open(partition_file, "r") as f:
        part = ""
        for l in f:
            l=l.rstrip("\n")
            if l[0] == "#":
                part = l.lstrip("#")
            else:
                read_id = l.split("\t")[0]
                partitions[part].add(read_id)

    return partitions