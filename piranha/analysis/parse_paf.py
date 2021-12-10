import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *

import os
import yaml

def group_hits(paf_file):

    hits = collections.defaultdict(list)
    with open(paf_file, "r") as f:
        for l in f:
            tokens = l.rstrip("\n").split("\t")
            hits[tokens[5]].append(tokens[0])
    
    ref_hits = {}
    not_mapped = []
    for hit in hits:
        if hit == '*':
            not_mapped = hits[hit]
        else:
            ref_hits[hit] = hits[hit]
    return ref_hits, not_mapped

def write_out_report_ref_reads(seq_index,ref_index,sample_composition,tax_outdir,hits,min_reads):
    refs_present = []
    with open(sample_composition,"w") as fw:
        fw.write("reference,num_reads,cns_build\n")
        for reference in hits:
            ref_hits = hits[reference]
            cns = False
            if len(ref_hits) >= min_reads:
                cns = True
                for enterovirus in ["Enterovirus","Echovirus","Coxsackievirus"]:
                    if reference.startswith(enterovirus):
                        cns = False

                if cns == True:
                    refs_present.append(reference)

                fw.write(f"{reference},{len(ref_hits)},{cns}\n")

                with open(os.path.join(tax_outdir,f"{reference}.reference.fasta"),"w") as fref:
                    ref_record = ref_index[reference]
                    fref.write(f">{ref_record.id}\n{ref_record.seq}\n")
                
                with open(os.path.join(tax_outdir,f"{reference}.reads.fastq"),"w") as freads:
                    for hit in ref_hits:
                        SeqIO.write(seq_index[hit],freads,"fastq")
            else:
                fw.write(f"{reference},{len(ref_hits)},{cns}\n")

def write_new_config(hits,config_out,barcode,config):
    barcode_config = config
    barcode_config["barcode"] = barcode
    barcode_config["taxa_present"] = list(hits.keys())
    with open(config_out, 'w') as fw:
        yaml.dump(barcode_config, fw) 

def parse_paf_file(paf_file,read_file,sample_composition,references_sequences,tax_outdir,config_out,barcode,config):

    hits, not_mapped = group_hits(paf_file)

    seq_index = SeqIO.index(read_file, "fastq")
    ref_index =  SeqIO.index(references_sequences,"fasta")

    write_out_report_ref_reads(seq_index,ref_index,sample_composition,tax_outdir,hits,config[KEY_MIN_READS])

    write_new_config(hits,config_out,barcode,config)


