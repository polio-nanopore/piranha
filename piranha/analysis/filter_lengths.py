import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *


def filter_reads_by_length(reads_in,reads_out,lengths_out,config):
    fastq_records = []
    length_list = []

    with open(reads_out,"w") as fw:
        for record in SeqIO.parse(reads_in,"fastq"):
            length = len(record)
            length_list.append(length)
            if length > int(config[KEY_MIN_READ_LENGTH]) and length < int(config[KEY_MAX_READ_LENGTH]):
                fastq_records.append(record)

        SeqIO.write(fastq_records,fw, "fastq")
    with open(lengths_out,"w") as fw:
        # this is in case we want to make a read length dist in the report
        for i in length_list:
            fw.write(f"{i}\n")
        fw.write("\n")
