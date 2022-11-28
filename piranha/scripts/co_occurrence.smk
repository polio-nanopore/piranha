import os
import collections
from Bio import SeqIO
import yaml
import json

from piranha.analysis.get_co_occurrence import *

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *



BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]


rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"co_occurrence_info.json")


rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        cns= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule minimap2_ref:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","mapped.ref.sam")
    shell:
        """
        minimap2 -ax map-ont --score-N=0 --secondary=no {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

rule minimap2_cns:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}_cns.minimap2.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","mapped.cns.sam")
    shell:
        """
        minimap2 -ax map-ont --score-N=0 --secondary=no {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

rule sam_to_seq:
    input:
        sam_ref = os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","mapped.ref.sam"),
        sam_cns = os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","mapped.cns.sam"),
        ref = rules.files.params.ref,
        cns = rules.files.params.cns
    params:
        reference = "{reference}"
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","pseudoaln.fasta")
    run:
        if "Sabin" in params.reference:
            shell("""
                gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam_ref:q} -o {output[0]:q} &> {log}
                """)
        else:
            shell("""
                gofasta sam toMultiAlign -r {input.cns:q} -s {input.sam_cns:q} -o {output[0]:q} &> {log}
                """)

rule get_co_occurrence:
    input:
        fasta = expand(os.path.join(config[KEY_TEMPDIR],"haplotype","{reference}","pseudoaln.fasta"),reference=REFERENCES),
        csv = os.path.join(config[KEY_TEMPDIR],"variants.csv")
    output:
        json = os.path.join(config[KEY_TEMPDIR],"co_occurrence_info.json")
    run:
        counter = {}
        with open(input.csv,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "Sabin" in row["reference"]:
                    variants = row["variants"]
                    ref_file = ""
                    for ref_file in input.fasta:
                        if row["reference"] in ref_file:

                            read_fasta_file = ref_file
                    if variants:
                        counter[row["reference"]] = get_combinations(variants,read_fasta_file,row["reference"],BARCODE,10)
        
        with open(output.json, "w") as fw:
            fw.write(json.dumps(counter))

# rule gather_reads_for_cooccurrence:
#     input:
#         json = os.path.join(config[KEY_TEMPDIR],"co_occurrence_info.json"),
#         reads = rules.files.params.reads
#     output:
