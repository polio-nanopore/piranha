#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *




def get_combinations(variants,read_fasta_file,reference,barcode,threshold):
    
    if variants:
        print(green(f"Co-occurrences for {reference} in {barcode}"))
        site_combinations = collections.Counter()
        variant_list = variants.split(";")
        sites = [int(i.split(":")[0]) for i in variant_list]

        c = 100
        for record in SeqIO.parse(read_fasta_file,KEY_FASTA):
            
            if not c:
                break

            bases = ";".join([f"{i}{record.seq[i-1]}" for i in sites])

            site_combinations[bases] +=1

            c -= 1

        sig_combinations = {}
        for i in site_combinations:
            
            pcentage  = site_combinations[i]
            if pcentage > threshold:
                
                print(i, site_combinations[i])
                sig_combinations[i] = pcentage
        
        return sig_combinations 

    else:
        return None 



# variants = "17:CT;160:GA"
# read_fasta_file = "/Users/s1680070/repositories/piranha/analysis_2022-07-15/barcode07/reference_analysis/Poliovirus3-Sabin_AY184221/pseudoaln.fasta"
# threshold = 5

# sig_combinations = get_combinations(variants,read_fasta_file,threshold)

# print(sig_combinations)
