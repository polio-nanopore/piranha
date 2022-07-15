#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv



def get_combinations(variants,read_fasta_file,threshold):
    
    if variants:

        site_combinations = collections.Counter()
        variant_list = variants.split(";")
        sites = [int(i.split(":")[0]) for i in variant_list]

        c = 1000
        for record in SeqIO.parse(read_fasta_file,"fasta"):
            
            if not c:
                break

            bases = ";".join([record.seq[i-1] for i in sites])
            print(bases)
            site_combinations[bases] +=1

            c -= 1

        sig_combinations = {}
        for i in site_combinations:
            print(i, site_combinations[i])
            pcentage  = site_combinations[i]
            if pcentage > threshold:
                sig_combinations[i] = pcentage
        
        return sig_combinations 

    else:
        return None 

# def get_co_occurring_reads():

variants = "17:CT;160:GA"
read_fasta_file = "/Users/s1680070/repositories/piranha/analysis_2022-07-15/barcode07/reference_analysis/Poliovirus3-Sabin_AY184221/pseudoaln.fasta"
threshold = 10

sig_combinations = get_combinations(variants,read_fasta_file,threshold)


