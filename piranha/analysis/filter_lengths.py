import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *


def filter_reads_by_length(reads_in,reads_out,config):
    fastq_records = []

    with open(reads_out,"w") as fw:
        for record in SeqIO.parse(reads_in,"fastq"):
            length = len(record)
            if length > int(config[KEY_MIN_READ_LENGTH]) and length < int(config[KEY_MAX_READ_LENGTH]):
                fastq_records.append(record)

        SeqIO.write(fastq_records,fw, "fastq")
    
