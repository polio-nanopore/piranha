import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *

import os
import yaml

def group_hits(paf_file):

    hits = collections.defaultdict(set)
    with open(paf_file, "r") as f:
        for l in f:
            tokens = l.rstrip("\n").split("\t")
            hits[tokens[5]].add(tokens[0])
    
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
            ref_hits = list(set(hits[reference]))
            cns = False
            if len(ref_hits) >= min_reads:
                cns = True
                for enterovirus in ["Enterovirus","Echovirus","Coxsackievirus"]:
                    if reference.startswith(enterovirus):
                        cns = False

                fw.write(f"{reference},{len(ref_hits)},{cns}\n")
                if cns == True:
                    refs_present.append(reference)

                    with open(os.path.join(tax_outdir,f"{reference}.fasta"),"w") as fref:
                        ref_record = ref_index[reference]
                        fref.write(f">{ref_record.id}\n{ref_record.seq}\n")
                    
                    with open(os.path.join(tax_outdir,f"{reference}.fastq"),"w") as freads:
                        for hit in ref_hits:
                            SeqIO.write(seq_index[hit],freads,"fastq")
            else:
                fw.write(f"{reference},{len(ref_hits)},{cns}\n")
    return refs_present

def write_new_config(hits,refs_present,config_out,barcode,tax_outdir,config):
    barcode_config = config
    barcode_config[KEY_BARCODE] = barcode
    barcode_config[KEY_TAXA_PRESENT] = refs_present
    barcode_config[KEY_HITS] = list(hits.keys())
    barcode_config[KEY_TAXA_OUTDIR] = tax_outdir

    with open(config_out, 'w') as fw:
        yaml.dump(barcode_config, fw) 

def parse_paf_file(paf_file,read_file,sample_composition,references_sequences,tax_outdir,config_out,barcode,config):

    hits, not_mapped = group_hits(paf_file)

    seq_index = SeqIO.index(read_file, "fastq")
    ref_index =  SeqIO.index(references_sequences,"fasta")

    refs_present = write_out_report_ref_reads(seq_index,ref_index,sample_composition,tax_outdir,hits,config[KEY_MIN_READS])

    write_new_config(hits,refs_present,config_out,barcode,tax_outdir,config)

    return refs_present


