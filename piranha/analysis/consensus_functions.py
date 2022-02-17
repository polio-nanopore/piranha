#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv
from itertools import groupby, count

from piranha.utils.config import *

def id_reference_cns(aln):
    ref_seq = False
    for record in SeqIO.parse(aln, "fasta"):
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
        bases = [query_seq[i],reference_seq[i]]
        if bases[0] != bases[1]:
            if bases[0] in non_amb and bases[1] in non_amb:
                snp = f"{i+1}:{bases[1]}{bases[0]}" # position-reference-query
                snps.append(snp)
            elif bases[0]=='-':
                insertions.append(i+1)
            elif bases[1]=='-':
                deletions.append(i+1)

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


def get_variation_pcent(ref,fasta):
    ref_record = SeqIO.read(ref,"fasta")
    
    #set up site counter, 0-based
    variant_sites = {}
    for i in range(len(ref_record)):
        variant_sites[i] = 0

    c = 0
    for record in SeqIO.parse(fasta,"fasta"):
        c +=1
        index = 0
        for pos in zip(str(record.seq),str(ref_record.seq)):
            
            if pos[0] != pos[1]:
                variant_sites[index]+=1
            index +=1

    variant_info = {}
    for site in variant_sites:
        variant_info[site] = collections.Counter()

    for record in SeqIO.parse(fasta,"fasta"):
        for site in variant_sites:
            variant = str(record.seq)[site]
            variant_info[site][variant]+=1

    x = []
    y = []
    info = []
    variation_info = []
    for site in variant_sites:
        
        pcent_variants = round(100*(variant_sites[site]/c), 1)
        x = site+1
        y = pcent_variants

        var_counts = variant_info[site]
        site_data = {"Position":x,"Percentage":y}
        for i in var_counts:
            site_data[f"{i} reads"] = var_counts[i]
        
        variation_info.append(site_data)

    return variation_info