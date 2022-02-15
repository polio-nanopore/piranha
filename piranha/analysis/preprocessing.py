import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os

def make_ref_display_name_map(references):
    ref_map = {}
    for record in SeqIO.parse(references,"fasta"):
        display_name = ""
        for item in str(record.description).split(" "):
            if item.startswith("display_name"):
                display_name = item.split("=")[1]
        ref_map[record.id] = display_name
    
    return ref_map

def group_hits(paf_file,ref_name_map):
    total_reads= 0
    hits = collections.defaultdict(set)
    with open(paf_file, "r") as f:
        for l in f:
            tokens = l.rstrip("\n").split("\t")
            hits[tokens[5]].add(tokens[0])
            total_reads +=1
    
    ref_hits = {}
    not_mapped = 0
    for hit in hits:
        if hit == '*':
            not_mapped += 1
        else:
            ref_hits[hit] = hits[hit]
    return ref_hits, not_mapped,total_reads

def write_out_report(ref_index,ref_map,csv_out,hits,unmapped,total_reads,barcode):

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

def write_out_hits(hits,outfile,min_reads):
    with open(outfile,"w") as fw:
        for hit in hits:
            reads = hits[hit]
            if len(reads) >= min_reads:
                for read in reads:
                    fw.write(f"{read},{hit}\n")

def parse_paf_file(paf_file,csv_out,hits_out,references_sequences,barcode,config):
    
    ref_name_map = make_ref_display_name_map(references_sequences)
    
    hits, unmapped,total_reads = group_hits(paf_file,ref_name_map)
    ref_index =  SeqIO.index(references_sequences,"fasta")
    write_out_report(ref_index,ref_name_map,csv_out,hits,unmapped,total_reads,barcode)

    write_out_hits(hits,hits_out,config[KEY_MIN_READS])


def diversity_report(input_files,csv_out,summary_out,barcodes_csv):
    summary_rows = {}
    with open(barcodes_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            summary_rows[row[KEY_BARCODE]] = {KEY_BARCODE: row[KEY_BARCODE],
                                            KEY_SAMPLE: row[KEY_SAMPLE]
                                            }
            for field in SAMPLE_SUMMARY_HEADER_FIELDS:
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
                    writer.writerow(row)
                    
    with open(summary_out,"w") as fw2:
        writer = csv.DictWriter(fw2, fieldnames=SAMPLE_SUMMARY_HEADER_FIELDS, lineterminator="\n")
        writer.writeheader()
        for barcode in summary_rows:
            row = summary_rows[barcode]

            writer.writerow(row)

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
            hit,ref = l.rstrip("\n").split(",")
            try:
                record = seq_index[hit]
                handle = handle_dict[ref]
                SeqIO.write(record,handle,"fastq")
            except:
                pass

    for ref in handle_dict:
        handle_dict[ref].close()



