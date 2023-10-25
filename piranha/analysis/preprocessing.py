#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os
import gzip
import io
import re

from piranha.utils.log_colours import green,cyan,red

def gather_filter_reads_by_length(dir_in,barcode,reads_out,config):

    if not os.path.exists(dir_in):
        os.mkdir(dir_in)

    fastq_records = []

    total_reads = 0
    with open(reads_out,"w") as fw:
        for r,d,f in os.walk(dir_in):
            for reads_in in f:
                if reads_in.endswith(".gz") or reads_in.endswith(".gzip"):
                    with gzip.open(os.path.join(dir_in,reads_in), "rt") as handle:
                        for record in SeqIO.parse(handle, "fastq"):
                            total_reads +=1
                            length = len(record)
                            if length > int(config[KEY_MIN_READ_LENGTH]) and length < int(config[KEY_MAX_READ_LENGTH]):
                                fastq_records.append(record)
                                

                elif reads_in.endswith(".fastq") or reads_in.endswith(".fq"):
                    for record in SeqIO.parse(os.path.join(dir_in,reads_in),"fastq"):
                        total_reads +=1
                        length = len(record)
                        if length > int(config[KEY_MIN_READ_LENGTH]) and length < int(config[KEY_MAX_READ_LENGTH]):
                            fastq_records.append(record)

        print(green(f"Total reads {barcode}:"),total_reads)
        print(green(f"Total passed reads {barcode}:"),len(fastq_records))
        SeqIO.write(fastq_records,fw, "fastq")

def make_ref_display_name_map(references):
    ref_map = {}
    for record in SeqIO.parse(references,KEY_FASTA):
        display_name = ""
        for item in str(record.description).split(" "):
            if item.startswith(KEY_DISPLAY_NAME):
                display_name = item.split("=")[1]
        ref_map[record.id] = display_name
    
    return ref_map

def make_display_name_to_reference_group_map(ref_map):
    ref_group_map = {}
    for key in ref_map:
        val_lower = ref_map[key].lower()
        if "wpv" in val_lower:
            num = re.findall(r'\d+', val_lower.split("wpv")[1])[0]
            ref_group_map[key] = "WPV%s" % num
        elif "nonpolio" in val_lower or "non-polio" in val_lower or "non_polio" in val_lower:
            ref_group_map[key] = "NonPolioEV"
        elif "polio" in val_lower and "wt" in val_lower:
            num = re.findall(r'\d+', val_lower.split("polio")[1])[0]
            ref_group_map[key] = "WPV%s" % num
        elif "polio" in val_lower:
            num = re.findall(r'\d+', val_lower.split("polio")[1])[0]
            ref_group_map[key] = "Sabin%s-related" % num
        elif "sabin" in val_lower:
            num = re.findall(r'\d+', val_lower.split("sabin")[1])[0]
            ref_group_map[key] = "Sabin%s-related" % num
        else:
            ref_group_map[key] = "NonPolioEV"
    assert len(ref_map) == len(ref_group_map)
    return ref_group_map


def parse_line(line):
    values = {}
    tokens = line.rstrip("\n").split("\t")
    values["read_name"], values["read_len"] = tokens[:2]

    values["read_hit_start"] = int(tokens[2])
    values["read_hit_end"] = int(tokens[3])
    values["direction"] = tokens[4]
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"],values["map_quality"] = tokens[5:12]
    
    values["ref_len"] = int(values["ref_len"])
    values["aln_block_len"] = int(values["aln_block_len"])

    return values

def add_to_hit_dict(hits, mapping,min_map_len,min_map_quality,unmapped):
    status,description = "",""
    if mapping["direction"] == "+":
        start = mapping["read_hit_start"]
        end = mapping["read_hit_end"]
    elif mapping["direction"] == "-":
        start = mapping["read_hit_end"]
        end = mapping["read_hit_start"]
    else:
        unmapped+=1
        status = "unmapped"

    if not status:
        if int(mapping["aln_block_len"]) > min_map_len:
            if int(mapping["map_quality"]) > min_map_quality:
                hits[mapping["ref_hit"]].add((mapping["read_name"],start,end,mapping["aln_block_len"]))
                status = "mapped"
            else:
                unmapped+=1
                status = "filtered"
                description = "mapping quality too low"
        else:
            unmapped+=1
            status = "filtered"
            description = "alignment block too short"

    return unmapped,status,description


def group_hits(paf_file,ref_name_map,min_aln_block,min_map_quality,mapping_filter_file):
    total_reads= 0
    ambiguous =0
    unmapped = 0
    hits = collections.defaultdict(set)
    
    mappings = []

    last_mapping = None
    with open(mapping_filter_file,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=["read_name","status","description"],lineterminator='\n')
        writer.writeheader()
        with open(paf_file, "r") as f:
            for l in f:
                filtered_reason = ""
                mapping = parse_line(l)
                total_reads +=1

                if last_mapping:
                    if mapping["read_name"] == last_mapping["read_name"]:
                        mappings.append(last_mapping)
                    else:
                        mappings.append(last_mapping)

                        if len(mappings) > 1:
                            ambiguous +=1
                            filtered_reason = "ambiguous mapping"
                            row = {"read_name":mapping["read_name"],"status":"unmapped","description":"ambiguous_mapping"}
                            writer.writerow(row)
                        else:
                            unmapped,status,description = add_to_hit_dict(hits, last_mapping,min_aln_block,min_map_quality,unmapped)
                            row = {"read_name":mapping["read_name"],"status":status,"description":description}
                            writer.writerow(row)

                        mappings = []
                        last_mapping = mapping

                else:
                    last_mapping = mapping

            unmapped,status,description = add_to_hit_dict(hits, last_mapping,min_aln_block,min_map_quality,unmapped)
            row = {"read_name":mapping["read_name"],"status":status,"description":description}
            writer.writerow(row)
            total_reads +=1
    
    if total_reads == 0:
        return {}, 0, 0, 0

    return hits, unmapped, ambiguous, total_reads

def write_out_report(ref_index,ref_map,csv_out,hits,unmapped,total_reads,barcode):

    with open(csv_out,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=SAMPLE_HIT_HEADER_FIELDS,lineterminator="\n")
        writer.writeheader()
        
        if total_reads == 0:
            pcent_unmapped = 100
        else:
            pcent_unmapped = round((100*(unmapped/total_reads)),2)

        unmapped_row = {
            KEY_BARCODE:barcode,
            KEY_REFERENCE:"unmapped",
            KEY_NUM_READS:unmapped,
            KEY_PERCENT: pcent_unmapped,
            KEY_REFERENCE_GROUP:"unmapped"}
        writer.writerow(unmapped_row)

        for reference in hits:
            hit_count = len(hits[reference])
            ref_group = ref_map[reference]
            mapped_row = {
                KEY_BARCODE:barcode,
                KEY_REFERENCE:reference,
                KEY_NUM_READS:hit_count,
                KEY_PERCENT: round((100*(hit_count/total_reads)),2),
                KEY_REFERENCE_GROUP:ref_group}
            writer.writerow(mapped_row)

def write_out_hits(hits,outfile):
    with open(outfile,"w") as fw:
        writer = csv.DictWriter(fw, lineterminator="\n",fieldnames=["read_name","hit","start","end","aln_block_len"])
        writer.writeheader()
        for hit in hits:
            hit_info = hits[hit]
            for read in hit_info:
                name,start,end,aln_len = read
                row = {"read_name":name,
                        "hit":hit,
                        "start":start,
                        "end":end,
                        "aln_block_len":aln_len}
                writer.writerow(row)

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def parse_paf_file(paf_file,
                    csv_out,
                    mapping_filter_out,
                    hits_out,
                    references_sequences,
                    barcode,
                    analysis_mode,
                    min_map_quality,
                    config):
    
    if is_non_zero_file(paf_file):
        
        permissive_ref_name_map = make_ref_display_name_map(references_sequences)
        ref_name_map = make_display_name_to_reference_group_map(permissive_ref_name_map)

        min_aln_block = config[KEY_MIN_ALN_BLOCK]

        ref_hits, unmapped, ambiguous, total_reads = group_hits(paf_file,ref_name_map,min_aln_block,min_map_quality,mapping_filter_out)
        print(f"Barcode: {barcode}")
        print(green("Unmapped:"),unmapped)
        print(green("Ambiguous mapping:"),ambiguous)

        ref_index =  SeqIO.index(references_sequences,KEY_FASTA)
        write_out_report(ref_index,ref_name_map,csv_out,ref_hits,unmapped,total_reads,barcode)

        write_out_hits(ref_hits,hits_out)
        
    else:
        print("No reads for",barcode)
        with open(csv_out,"w") as fw:
            writer = csv.DictWriter(fw, fieldnames=SAMPLE_HIT_HEADER_FIELDS,lineterminator="\n")
            writer.writeheader()

        with open(hits_out,"w") as fw:
            writer = csv.DictWriter(fw, lineterminator="\n",fieldnames=["read_name","hit","start","end","aln_block_len"])
            writer.writeheader()

        with open(mapping_filter_out,"w") as fw:
            writer = csv.DictWriter(fw, lineterminator="\n",fieldnames=["read_name","status","description"])
            writer.writeheader()



def diversity_report(input_files,csv_out,summary_out,ref_file,config):
    barcodes_csv= config[KEY_BARCODES_CSV]
    min_reads = config[KEY_MIN_READS]
    min_pcent = config[KEY_MIN_PCENT]
    summary_rows = {}
    refs_out = collections.defaultdict(list)

    SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1

    if config[KEY_ANALYSIS_MODE] == VALUE_ANALYSIS_MODE_WG:
        SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG
    

    with open(barcodes_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            summary_rows[row[KEY_BARCODE]] = {KEY_BARCODE: row[KEY_BARCODE],
                                            KEY_SAMPLE: row[KEY_SAMPLE]
                                            }
            for field in SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS:
                if field not in summary_rows[row[KEY_BARCODE]]: # the rest are ref counters
                    summary_rows[row[KEY_BARCODE]][field] = 0

    with open(csv_out,"w") as fw:
        writer = csv.DictWriter(fw,fieldnames=SAMPLE_HIT_HEADER_FIELDS,lineterminator="\n")
        writer.writeheader()

        for report_file in input_files:
            barcode = ""
            with open(report_file, "r") as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    if not barcode:
                        barcode = row[KEY_BARCODE]

                    summary_rows[barcode][row[KEY_REFERENCE_GROUP]] += int(row[KEY_NUM_READS])

                    if int(row[KEY_NUM_READS]) >= min_reads and float(row[KEY_PERCENT]) >= min_pcent:
                        refs_out[barcode].append(row[KEY_REFERENCE])
                        writer.writerow(row)


    with open(summary_out,"w") as fw2:
        writer = csv.DictWriter(fw2, fieldnames=SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS, lineterminator="\n")
        writer.writeheader()
        for barcode in summary_rows:
            row = summary_rows[barcode]
            writer.writerow(row)

    config[KEY_BARCODES] = []
    for barcode in refs_out:
        refs = refs_out[barcode]
        
        refs = [i for i in refs if i != "unmapped"]
        config[barcode]=refs
        if refs:
            config[KEY_BARCODES].append(barcode)
    
    return config

def check_which_refs_to_write(input_csv,min_reads,min_pcent):
    to_write = set()
    with open(input_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row[KEY_NUM_READS]) >= min_reads and float(row[KEY_PERCENT]) >= min_pcent:
                if row[KEY_REFERENCE] != "unmapped":
                    to_write.add(row[KEY_REFERENCE])
                    print(f"{row[KEY_REFERENCE_GROUP]}\t{row[KEY_NUM_READS]} reads\t{row[KEY_PERCENT]}% of sample")
    return list(to_write)

def write_out_fastqs(input_csv,input_hits,input_fastq,outdir,primer_length,config):
    seq_index = SeqIO.index(input_fastq, "fastq")
    
    to_write = check_which_refs_to_write(input_csv,config[KEY_MIN_READS],config[KEY_MIN_PCENT])
    handle_dict = {}
    for ref in to_write:
        handle_dict[ref] = open(os.path.join(outdir,f"{ref}.fastq"),"w")
    
    not_written = collections.Counter()

    with open(input_hits,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                read_name = row["read_name"]
                record = seq_index[read_name]

                hit = row["hit"]
                handle = handle_dict[hit]

                trimmed_record = record[primer_length:-primer_length]

                SeqIO.write(trimmed_record,handle,"fastq")
            except:
                not_written[hit]+=1

    # if not_written:
    #     print("Read hits not written because below threshold:")
    #     for i in not_written:
    #         print(i, not_written[i])

    for ref in handle_dict:
        handle_dict[ref].close()

    return to_write

def write_out_ref_fasta(to_write,ref_file,outdir):
    ref_index = SeqIO.index(ref_file, KEY_FASTA)

    for ref in to_write:
        with open(os.path.join(outdir,f"{ref}.reference.fasta"),"w") as fw:
            record = ref_index[ref]
            fw.write(f">{record.description}\n{record.seq}\n")
    
