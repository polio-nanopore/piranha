#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv
from itertools import groupby, count
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

# def adjust_position(cs,primer_length):

#     # new cigar string
#     c = Cigar(cs)
#     items = list(c.items())
#     start_mask = items[0][0]
#     end_mask = items[-1][0]

#     amended_start_mask = start_mask + primer_length
#     amended_end_mask = end_mask + primer_length

#     new_cigar = c.mask_left(amended_start_mask).mask_right(amended_end_mask)

#     return new_cigar.cigar

# def soft_mask_primer_sites(input_sam, output_sam, primer_length):
#     with open(output_sam, "w") as fw:
#         with open(input_sam,"r") as f:
#             for l in f:
#                 l=l.rstrip("\n")

#                 if not l.startswith("@"):
#                     tokens = l.split("\t")
#                     new_tokens = tokens
#                     new_tokens[5] = adjust_position(tokens[5],primer_length)
#                     new_l = "\t".join(new_tokens)
#                     fw.write(f"{new_l}\n")
#                 else:
#                     fw.write(f"{l}\n")


#Get the reference bases at each position and stick it in a dictionary
def ref_dict_maker(ref_fasta):
    ref_dict = {}
    ref_fasta = pysam.FastaFile(ref_fasta)
    for idx, base in enumerate(ref_fasta.fetch(ref_fasta.references[0])):
        ref_dict[idx] = base

    return ref_dict

#function to calculate the percentage of non ref bases at a given position
def non_ref_prcnt_calc(pos,pileup_dict,ref_dict):
    non_ref_prcnt = 0
    ref_count = 0
    total = 0
    if ref_dict[pos] != "N":
        ref_base = ref_dict[pos] + " reads"
        for key in pileup_dict:
            if key != "Position":
                total += pileup_dict[key]
        ref_count = pileup_dict[ref_base]
        non_ref_prcnt = round((100 - ((ref_count / total) * 100)), 2)

    return non_ref_prcnt

#Use mpileup to get bases per read at each postion, then calculate % vs ref for each
def pileupper(bamfile,ref_dict,base_q=13):
    bamfile = pysam.AlignmentFile(bamfile, "rb" )
    variation_info = []
    for pileupcolumn in bamfile.pileup(bamfile.references[0], min_base_quality=base_q):
        pileup_dict = {}

        A_counter, G_counter, C_counter, T_counter, del_counter = 0,0,0,0,0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                counts_list = []
                #print(pileupread.alignment.query_position)
                #print(pileupread.alignment.query_sequence[pileupread.query_position])
                if pileupread.alignment.query_sequence[pileupread.query_position] == "A":
                    A_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "G":
                    G_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "C":
                    C_counter += 1
                elif pileupread.alignment.query_sequence[pileupread.query_position] == "T":
                    T_counter += 1
            else:
                del_counter += 1
        pileup_dict["Position"] = pileupcolumn.pos + 1
        pileup_dict["A reads"] = A_counter
        pileup_dict["C reads"] = C_counter
        pileup_dict["T reads"] = T_counter
        pileup_dict["G reads"] = G_counter
        pileup_dict["- reads"] = del_counter
        pileup_dict["Percentage"] = non_ref_prcnt_calc(pileupcolumn.pos,pileup_dict,ref_dict)
        pileup_dict["ref_base"] = ref_dict[pileupcolumn.pos]
        variation_info.append(pileup_dict)

    return (variation_info)
