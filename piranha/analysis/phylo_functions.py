#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

from piranha.utils.log_colours import green,cyan,red



def get_seqs_and_clusters(sample_seqs,supplementary_sequences,reference_sequences,phylo_outdir,config):
    seq_clusters = collections.defaultdict(list)
    for record in SeqIO.parse(sample_seqs,"fasta"):
        for ref_group in config[KEY_REFERENCES_FOR_CNS]:
            if ref_group in record.id:
                new_record = record
                new_record.id = "|".join(new_record.id.split("|")[:2])
                seq_clusters[ref_group].append(record)

    print(green("Reference groups for phylo pipeline:"))
    for i in seq_clusters:
        print(f"- {i}")
    
    if supplementary_sequences:
        for record in SeqIO.parse(supplementary_sequences,"fasta"):
            for ref_group in seq_clusters:
                if ref_group in record.description:
                    seq_clusters[ref_group].append(record)
    
    for record in SeqIO.parse(reference_sequences, "fasta"):
        for ref_group in seq_clusters:
            if ref_group in record.description:
                seq_clusters[ref_group].append(record)

    for i in seq_clusters:
        print(i, len(seq_clusters[i]))

        with open(os.path.join(phylo_outdir, f"{i}.fasta"),"w") as fw:
            SeqIO.write(seq_clusters[i], fw, "fasta")

    return list(seq_clusters.keys())

