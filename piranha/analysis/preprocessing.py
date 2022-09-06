#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

def make_ref_display_name_map(references):
    ref_map = {}
    for record in SeqIO.parse(references,KEY_FASTA):
        display_name = ""
        for item in str(record.description).split(" "):
            if item.startswith("display_name"):
                display_name = item.split("=")[1]
        ref_map[record.id] = display_name
    
    return ref_map

def group_hits(paf_file,padding,ref_name_map):
    total_reads= 0
    hits = collections.defaultdict(set)
    with open(paf_file, "r") as f:
        for l in f:
            tokens = l.rstrip("\n").split("\t")
            read_hit_start = int(tokens[2])+padding
            read_hit_end = int(tokens[3])-padding
            if tokens[4] == "+":
                hits[tokens[5]].add((tokens[0],read_hit_start,read_hit_end))
            else:
                hits[tokens[5]].add((tokens[0],read_hit_end,read_hit_start))
            total_reads +=1
    
    ref_hits = {}
    not_mapped = 0
    for hit in hits:
        if hit == '*':
            not_mapped += 1
        else:
            ref_hits[hit] = hits[hit]
    if total_reads == 0:
        return {}, 0, 0
    return ref_hits, not_mapped,total_reads

def write_out_report(ref_index,ref_map,csv_out,hits,unmapped,total_reads,barcode):

    if total_reads == 0:
        pcent_unmapped = 100
    else:
        pcent_unmapped = round((100*(unmapped/total_reads)),2)
    
    with open(csv_out,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=SAMPLE_HIT_HEADER_FIELDS,lineterminator="\n")
        writer.writeheader()
        
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
        for hit in hits:
            reads = hits[hit]
            for read,start,end in reads:
                fw.write(f"{read},{hit},{start},{end}\n")

def parse_paf_file(paf_file,csv_out,hits_out,references_sequences,barcode,analysis_mode,config):
    
    # this is how much of the mapped coords to mask from the sequences
    padding = 0
    if analysis_mode == VALUE_ANALYSIS_MODE_WG_2TILE:
        padding = 30

    ref_name_map = make_ref_display_name_map(references_sequences)
    
    hits, unmapped,total_reads = group_hits(paf_file,padding,ref_name_map)
    ref_index =  SeqIO.index(references_sequences,KEY_FASTA)
    write_out_report(ref_index,ref_name_map,csv_out,hits,unmapped,total_reads,barcode)

    write_out_hits(hits,hits_out)


def diversity_report(input_files,csv_out,summary_out,config):
    barcodes_csv= config[KEY_BARCODES_CSV]
    min_reads = config[KEY_MIN_READS]
    min_pcent = config[KEY_MIN_PCENT]
    summary_rows = {}
    refs_out = collections.defaultdict(list)

    SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS = SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1

    if config[KEY_ANALYSIS_MODE] == VALUE_ANALYSIS_MODE_WG_2TILE:
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
            with open(report_file, "r") as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    summary_rows[row[KEY_BARCODE]][row[KEY_REFERENCE_GROUP]] += int(row[KEY_NUM_READS])
                    if int(row[KEY_NUM_READS]) >= min_reads and float(row[KEY_PERCENT]) >= min_pcent:
                        refs_out[row[KEY_BARCODE]].append(row[KEY_REFERENCE])
                        writer.writerow(row)
                    
    with open(summary_out,"w") as fw2:
        writer = csv.DictWriter(fw2, fieldnames=SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS, lineterminator="\n")
        writer.writeheader()
        for barcode in summary_rows:
            row = summary_rows[barcode]
            writer.writerow(row)

    config[KEY_BARCODES] = []
    for barcode in refs_out:
        refs = list(set(refs_out[barcode]))
        config[barcode]=refs
        config[KEY_BARCODES].append(barcode)
    
    return config

def check_which_refs_to_write(input_csv,min_reads,min_pcent):
    to_write = set()
    with open(input_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row[KEY_NUM_READS]) >= min_reads and float(row[KEY_PERCENT]) >= min_pcent:
                to_write.add(row[KEY_REFERENCE])
                print(f"{row[KEY_REFERENCE_GROUP]}\t{row[KEY_NUM_READS]} reads\t{row[KEY_PERCENT]}% of sample")
    return list(to_write)

def write_out_fastqs(input_csv,input_hits,input_fastq,outdir,config):
    seq_index = SeqIO.index(input_fastq, "fastq")
    
    to_write = check_which_refs_to_write(input_csv,config[KEY_MIN_READS],config[KEY_MIN_PCENT])

    handle_dict = {}
    for ref in to_write:
        handle_dict[ref] = open(os.path.join(outdir,f"{ref}.fastq"),"w")
    
    not_written = collections.Counter()
    with open(input_hits,"r") as f:
        for l in f:
            try:
                hit,ref,start,end = l.rstrip("\n").split(",")
                record = seq_index[hit]
                
                # trimmed = record[:int(start)].lower()+ record[int(start):int(end)] + record[:int(start)].lower()
                handle = handle_dict[ref]

                SeqIO.write(record,handle,"fastq")
            except:
                pass


    for ref in handle_dict:
        handle_dict[ref].close()
    
    return to_write

def write_out_ref_fasta(to_write,ref_file,outdir):
    ref_index = SeqIO.index(ref_file, KEY_FASTA)

    for ref in to_write:
        with open(os.path.join(outdir,f"{ref}.reference.fasta"),"w") as fw:
            record = ref_index[ref]
            fw.write(f">{record.description}\n{record.seq}\n")
    


