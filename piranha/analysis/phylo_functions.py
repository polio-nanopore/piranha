#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

from piranha.utils.log_colours import green,cyan,red


def get_seqs_and_clusters(sample_seqs,supplementary_sequences,reference_sequences,outgroup_sequences,barcodes_csv,phylo_outdir,config):
    

    seq_metadata = collections.defaultdict(dict)
    seq_clusters = collections.defaultdict(list)
    for record in SeqIO.parse(sample_seqs,"fasta"):
        for ref_group in config[KEY_REFERENCES_FOR_CNS]:
            if ref_group in record.id:
                new_record = record
                new_header = "|".join(new_record.id.split("|")[:2])
                new_record.description = new_header
                new_record.id = new_header
                seq_clusters[ref_group].append(new_record)

                seq_metadata[new_record.id]["sample"] = new_record.id
                seq_metadata[new_record.id]["query_boolean"] = "True"
                seq_metadata[new_record.id]["reference_group"] = ref_group


    print(green("Reference groups for phylo pipeline:"))
    for i in seq_clusters:
        print(f"- {i}")
    
    if supplementary_sequences:
        for record in SeqIO.parse(supplementary_sequences,"fasta"):
            for ref_group in seq_clusters:
                if ref_group in record.description:
                    seq_clusters[ref_group].append(record)

                    seq_metadata[record.id]["sample"] = record.id
                    seq_metadata[record.id]["query_boolean"] = "False"
                    seq_metadata[record.id]["reference_group"] = ref_group
    
    for record in SeqIO.parse(reference_sequences, "fasta"):
        for ref_group in seq_clusters:
            if ref_group in record.description:
                seq_clusters[ref_group].append(record)

                seq_metadata[record.id]["sample"] = record.id
                seq_metadata[record.id]["query_boolean"] = "False"
                seq_metadata[record.id]["reference_group"] = ref_group
    
    for record in SeqIO.parse(outgroup_sequences, "fasta"):
        for ref_group in seq_clusters:
            if ref_group in record.description:
                new_record = record
                new_record.id = "outgroup"
                new_record.description = "outgroup"
                seq_clusters[ref_group].append(new_record)

    header = ["sample","query_boolean","reference_group"]
    with open(barcodes_csv, "r") as f:
        reader = csv.DictReader(f)
        reader_header = reader.fieldnames
        for col in reader_header:
            if col in config[KEY_PHYLO_METADATA_COLUMNS]:
                header.append(col)
        for row in reader:
            for col in header:
                seq_metadata[row[KEY_SAMPLE]][col] = row[col]

    with open(os.path.join(phylo_outdir, f"annotations.csv"), "w") as fw0:
        
        writer0 = csv.DictWriter(fw0,fieldnames=header, lineterminator='\n')
        writer0.writeheader()

        for i in seq_clusters:
            print(i, len(seq_clusters[i]))

            with open(os.path.join(phylo_outdir, f"{i}.annotations.csv"), "w") as fw:
                
                writer = csv.DictWriter(fw,fieldnames=header, lineterminator='\n')
                writer.writeheader()
                for record in seq_metadata:
                    if seq_metadata[record]["reference_group"] == i:
                        row = seq_metadata[record]
                        for col in header:
                            if col not in row:
                                row[col] = ""
                            new_data = row[col].replace("'","").replace(";","").replace("(","").replace(")","")
                            row[col] = new_data
                        writer.writerow(row)
                        writer0.writerow(row)

            with open(os.path.join(phylo_outdir, f"{i}.fasta"),"w") as fw:
                SeqIO.write(seq_clusters[i], fw, "fasta")

    return list(seq_clusters.keys())

