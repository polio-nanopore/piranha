#!/usr/bin/env python3

from Bio import SeqIO
from Bio import AlignIO
import sys

from piranha.analysis.consensus_functions import find_variants,id_reference_cns

def trim_trailing_gaps(alignment):
    for col in range(alignment.get_alignment_length()):
        if not "-" in alignment[:, col]:
            start_position = col
            break

    for col in range(alignment.get_alignment_length()):
        end_index = col+1
        if not "-" in alignment[:, -end_index]:
            end_position = col
            break

    print(f"\nTrimming trailing gaps in alignment.\nAlignment now from {start_position} to {alignment.get_alignment_length()-end_position}.\n")   

    if end_position == 0:
        return alignment[:,start_position:]
    else:
        return alignment[:,start_position:-end_position]

def remove_gaps(aln):

    untrimmed_alignment = AlignIO.read(aln, "fasta")
    trimmed = trim_trailing_gaps(untrimmed_alignment)
    
    cleaned_variants = []

    print(f"Reading in {aln}.\n\nGaps found:")
    cns_string = ""
    for i in range(len(trimmed[0])):
        col = trimmed[:,i]
        if len(set(col))>1 and '-' in col:
            print(f"Position {i+1}:\tReference:\t{col[0]}\tConsensus:\t{col[1]}")
            if col[0] == '-':
                pass
            else:
                cns_string += 'N'
        else:
            cns_string+= col[1]

    return trimmed[1].id, cns_string


def clean_cns_gaps(step, sample, aln_file, outfile):
    sample = sample.replace(" ","_")
    round_string = f" note={step}"

    cns_id, new_consensus = remove_gaps(aln_file)

    #the rule is to replace a gap in the query with 'N' and to force delete a base that causes a gap in the reference
    with open(outfile, "w") as fw:

        fw.write(f">{sample} accession={cns_id.split(':')[0]}{round_string} length={len(new_consensus)}\n{new_consensus.upper()}\n")

def clean_cns_mask(in_cns,to_mask):
    new_cns = ""
    masking = 0
    skipping = 0
    for i in range(len(in_cns)):
        if masking>0:
            new_cns += "N"
            masking -= 1
        elif skipping >0:
            skipping -= 1
        else:
            base = in_cns[i]
            if i in to_mask:
                variant = to_mask[i]
                var = variant.split(":")[1]
                if var.startswith("del"):
                    skipping = 0
                    masking = int(var[3:])
                    new_cns += "N"
                    masking -= 1
                elif var.startswith("ins"):
                    masking = 0
                    skipping = int(var[3:])
                    skipping -= 1
                else:
                    new_cns += "N"
            else:
                new_cns += base
    return new_cns

def clean_medaka_cns(sample, aln_file, outfile):
    sample = sample.replace(" ","_")

    ref_seq,cns_seq = id_reference_cns(aln_file)
    variants = find_variants(ref_seq,cns_seq)
    
    to_mask = {}
    window_dict = {}
    for variant in variants:

        pos,var = variant.split(":")
        if "del" in var or "ins" in var:
            length_indel = int(var[3:])
            if length_indel %3 != 0:
                index = int(pos)-1
                to_mask[index] = variant
        
        window_dict[variant] = 0
        for window_variant in window_dict:
            window_pos = int(window_variant.split(":")[0])
            window_range = range(window_pos-10,window_pos+10)
            if pos in window_range:
                window_dict[window_variant] +=1
        
    for variant in window_dict:
        if window_dict[variant] > 8:
            index = int(variant.split(":")[0])-1
            to_mask[index] = variant

    new_cns = clean_cns_mask(cns_seq,to_mask)

    with open(outfile, "w") as fw:

        fw.write(f">{sample} round=medaka_cleaned length={len(new_cns)}\n{new_cns.upper()}\n")

    return to_mask





