#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv
from itertools import groupby, count
import pysam
import pandas as pd
import json
# from cigar import Cigar

from piranha.utils.config import *
 
def id_reference_cns(aln): 
    ref_seq = False
    for record in SeqIO.parse(aln, KEY_FASTA):
        if not ref_seq:
            ref_seq = str(record.seq).upper()
        else:
            cns_seq = str(record.seq).upper()
    
    return ref_seq,cns_seq

def merge_indels(indel_list,prefix):
    if indel_list:
        groups = groupby(indel_list, key=lambda item, c=count():item-next(c))
        tmp = [list(g) for k, g in groups]
        merged_indels = []
        for i in tmp:
            indel = f"{i[0]}:{prefix}{len(i)}"
            merged_indels.append(indel)
        return merged_indels
    
    return indel_list

def find_variants(reference_seq,query_seq):

    non_amb = ["A","T","G","C"]
    
    variants = []

    snps =[]
    insertions = []
    deletions = []

    for i in range(len(query_seq)):
        # bases[0] == query seq
        # bases[1] == ref seq
        bases = [query_seq[i],reference_seq[i]]
        if bases[0] != bases[1]:
            
            if bases[0] in non_amb and bases[1] in non_amb:
                #if neither of them are ambiguous
                snp = f"{i+1}:{bases[1]}{bases[0]}" # position-reference-query
                snps.append(snp)
            elif bases[0]=='-':
                #if there's a gap in the query, means a deletion
                deletions.append(i+1)
            elif bases[1]=='-':
                #if there's a gap in the ref, means an insertion
                insertions.append(i+1)

    insertions = merge_indels(insertions,"ins")
    deletions = merge_indels(deletions,"del")

    for var_list in [snps,insertions,deletions]:
        for var in var_list:
            variants.append(var)

    variants = sorted(variants, key = lambda x : int(x.split(":")[0]))

    return variants

def find_ambiguity_pcent(query_seq):
    ambiguity_count = 0
    for i in query_seq:
        if i not in ["A","T","G","C","-"]:
            ambiguity_count +=1
    ambiguity_pcent = round(100*(ambiguity_count/len(query_seq)),2)
    return ambiguity_pcent

def parse_variants(alignment,out_report,barcode,reference):

    ref_seq,cns_seq = id_reference_cns(alignment)
    variants = find_variants(ref_seq,cns_seq)
    var_string = ";".join(variants)
    with open(out_report,"w") as fw:
        fw.write(f"{barcode},{reference},{len(variants)},{var_string}\n")


def join_variant_files(header_fields,in_files,output):
    with open(output,"w") as fw:
        header = ",".join(header_fields) + "\n"
        fw.write(header)
        for in_file in in_files:
            with open(in_file, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    fw.write(f"{l}\n")

