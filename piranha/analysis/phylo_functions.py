#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

from piranha.utils.log_colours import green,cyan,red


def get_seqs_and_clusters(sample_seqs,supplementary_sequences,reference_sequences,outgroup_sequences,barcodes_csv,supplementary_metadata,phylo_outdir,config):
    
    seq_metadata = collections.defaultdict(dict)
    seq_clusters = collections.defaultdict(list)

    header = VALUE_PHYLO_HEADER
    written = {}
    for record in SeqIO.parse(sample_seqs,KEY_FASTA):
        for ref_group in config[KEY_REFERENCES_FOR_CNS]:

            """
            record header is:
            >SAMPLE|REFERENCE_GROUP|CNS_ID|EPID|DATE barcode=barcode01 variant_count=8 variants=17:CT;161:CT;427:GA;497:AC;507:CT;772:AG;822:CT;870:CA 

            if "all_metadata" then everything else gets added to the description
            """
            
            fields = record.description.split(" ")
            
            record_sample,reference_group,cns_id,epid,sample_date = record.id.split("|")

            description_dict = {}
            for field in fields[1:]:
                key,value = field.split("=")
                description_dict[key] = value

            if ref_group == reference_group:
                
                new_record = record
                new_record.description = ""
                new_record.id = record.id

                barcode = description_dict[KEY_BARCODE]
                name = new_record.id


                seq_metadata[name][KEY_NAME] = name
                seq_metadata[name][KEY_SAMPLE] = record_sample
                seq_metadata[name][KEY_BARCODE] = barcode
                seq_metadata[name][KEY_SOURCE] = "Sample"
                seq_metadata[name][KEY_REFERENCE_GROUP] = ref_group

                var_count = description_dict[KEY_VARIANT_COUNT]

                call = ref_group
                if ref_group.startswith("Sabin"):
                    # configured number of mutations in sabin for the call threshold of VDPV
                    call_threshold = CALL_THRESHOLD_DICT[ref_group]
                    if int(var_count) > call_threshold:
                        call = "VDPV"
                    elif var_count == 0:
                        call = "Sabin"
                    else:
                        call = "Sabin-like"

                seq_metadata[name][KEY_CALL] = call
                seq_clusters[ref_group].append(new_record)
                continue


    print(green("Reference groups for phylo pipeline:"))
    for i in seq_clusters:
        print(f"- {i}")
    
    if supplementary_sequences:
        for record in SeqIO.parse(supplementary_sequences,KEY_FASTA):
            if record:
                for ref_group in seq_clusters:
                    if ref_group in record.description:
                        

                        seq_metadata[record.id][KEY_NAME] = record.id
                        seq_metadata[record.id][KEY_SAMPLE] = record.id
                        seq_metadata[record.id][KEY_SOURCE] = "Background"
                        seq_metadata[record.id][KEY_REFERENCE_GROUP] = ref_group
                        seq_metadata[record.id][KEY_CALL] = ref_group

                        new_record = record
                        new_record.description = ""
                        seq_clusters[ref_group].append(new_record)
    
    if supplementary_metadata:
        with open(supplementary_metadata, "r") as f:
            reader = csv.DictReader(f)
            for col in reader.fieldnames:
                if col in config[KEY_SUPPLEMENTARY_METADATA_COLUMNS] and col not in header:
                    header.append(col)

            for row in reader:
                if row[config[KEY_SUPPLEMENTARY_METADATA_ID_COLUMN]] in seq_metadata:
                    sample = row[config[KEY_SUPPLEMENTARY_METADATA_ID_COLUMN]]
                    for col in config[KEY_SUPPLEMENTARY_METADATA_COLUMNS]:
                        if col in row:
                            seq_metadata[sample][col] = row[col]

    for record in SeqIO.parse(reference_sequences, KEY_FASTA):
        if record:
            for ref_group in seq_clusters:
                if ref_group in record.description:
                    seq_clusters[ref_group].append(record)

                    seq_metadata[record.id][KEY_NAME] = record.id
                    seq_metadata[record.id][KEY_SAMPLE] = record.id
                    
                    seq_metadata[record.id][KEY_REFERENCE_GROUP] = ref_group
                    seq_metadata[record.id][KEY_CALL] = ref_group

                    if "Sabin" in record.description:
                        seq_metadata[record.id][KEY_SOURCE] = "Sabin"
                    else:
                        seq_metadata[record.id][KEY_SOURCE] = "Reference"
    
    for record in SeqIO.parse(outgroup_sequences, KEY_FASTA):
        for ref_group in seq_clusters:
            if ref_group in record.description:
                new_record = record
                new_record.description = ""
                new_record.id = "outgroup"
                
                seq_clusters[ref_group].append(new_record)

    with open(barcodes_csv, "r") as f:
        reader = csv.DictReader(f)
        reader_header = reader.fieldnames
        for col in reader_header:
            if col in config[KEY_PHYLO_METADATA_COLUMNS] and col not in header:
                header.append(col)
        
        for row in reader:
            for k in seq_metadata:
                if k.split("|")[0] ==row[KEY_SAMPLE]:
                    print(k, row[KEY_SAMPLE])
                    for col in header:
                        if col in reader_header and col not in [KEY_SAMPLE,KEY_BARCODE,KEY_SOURCE,KEY_NAME]:
                            seq_metadata[k][col] = row[col]
                    print(seq_metadata[k])

    with open(os.path.join(phylo_outdir, f"annotations.csv"), "w") as fw0:
        
        writer0 = csv.DictWriter(fw0,fieldnames=header, lineterminator='\n')
        writer0.writeheader()

        for i in seq_clusters:
            print(i, len(seq_clusters[i]))

            with open(os.path.join(phylo_outdir, f"{i}.annotations.csv"), "w") as fw:
                
                writer = csv.DictWriter(fw,fieldnames=header, lineterminator='\n')
                writer.writeheader()
                for record in seq_metadata:
                    if seq_metadata[record][KEY_REFERENCE_GROUP] == i:
                        row = seq_metadata[record]
                        for col in header:
                            if col not in row:
                                row[col] = ""
                            new_data = str(row[col]).replace("'","").replace(";","").replace("(","").replace(")","")
                            row[col] = new_data
                        writer.writerow(row)
                        writer0.writerow(row)

            with open(os.path.join(phylo_outdir, f"{i}.fasta"),"w") as fw:
                records = seq_clusters[i]
                record_dict = {}
                for record in records:
                    record_dict[record.id] = record
                unique_records = [r for r in record_dict.values()]
                SeqIO.write(unique_records, fw, KEY_FASTA)

    tree_annotations = config[KEY_TREE_ANNOTATIONS]
    for i in header:
        if i not in [KEY_SAMPLE,KEY_BARCODE,KEY_SOURCE,KEY_NAME]:
            tree_annotations+= f"{i} "


    return list(seq_clusters.keys()),tree_annotations

def update_local_database(sample_sequences,detailed_csv,new_db_seqs,new_db_metadata,config):
    
    record_ids = {}
    with open(new_db_seqs,"w") as fw:
        countnew = 0

        for record in SeqIO.parse(sample_sequences, KEY_FASTA):
            new_record = record
            desc_list = new_record.description.split(" ")
            write_record = True

            for i in desc_list:
                if i.startswith("variant_count"):
                    count = i.split("=")[1]
                elif i.startswith(VALUE_REFERENCE_GROUP_FIELD):
                    match_field = i.split("=")[1]

            if "Sabin" in match_field:
                if int(count) < config[KEY_LOCAL_DATABASE_THRESHOLD]:
                    write_record = False
            elif match_field == "PositiveControl":
                write_record = False
            
            if write_record:
                new_desc_list = [i for i in desc_list if not i.startswith("barcode=")]
                new_record.description = " ".join(new_desc_list)
                
                SeqIO.write(new_record, fw, KEY_FASTA)
                countnew+=1
                sample = record.id.split("|")[0]
                record_ids[record.id] = sample

    with open(new_db_metadata,"w") as fw:
        with open(detailed_csv,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            header.append(config[KEY_SUPPLEMENTARY_METADATA_ID_COLUMN])

            writer = csv.DictWriter(fw, fieldnames=header, lineterminator="\n")
            writer.writeheader()
            sample_data = {}
            for row in reader:
                sample = row[KEY_SAMPLE]
                sample_data[sample] = row

            for record_id in record_ids:
                sample = record_ids[record_id]
                row = sample_data[sample]
                row[config[KEY_SUPPLEMENTARY_METADATA_ID_COLUMN]] = record_id
                writer.writerow(row)

    print(green(f"Local database updated with ")+ f"{countnew}"+ green(" newly generated records."))
