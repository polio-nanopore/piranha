#!/usr/bin/env python3

from Bio import SeqIO
from Bio import AlignIO
import sys


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


def clean_racon_gaps(round_count, tax_name, aln_file, outfile):
    round_string = f" round_name={round_count}"

    cns_id, new_consensus = remove_gaps(aln_file)


    #the rule is to replace a gap in the query with 'N' and to force delete a base that causes a gap in the reference
    with open(outfile, "w") as fw:

        fw.write(f">{tax_name} accession={cns_id.split(':')[0]}{round_string} length={len(new_consensus)}\n{new_consensus.upper()}\n")
