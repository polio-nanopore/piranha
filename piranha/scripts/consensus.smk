import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

BARCODE = config[KEY_BARCODE]
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta"),
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus.fasta"), reference=REFERENCES)

rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule minimap2:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam")
    shell:
        """
        minimap2 -ax map-ont {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

rule sort_index:
    input:
        sam = rules.minimap2.output.sam
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.sorted.bam"),
        index = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.sorted.bam.bai")
    shell:
        """
        samtools view -bS -F 4 {input.sam:q} | samtools sort -o {output[0]:q} &&
        samtools index {output.bam:q} {output.index:q}
        """

rule medaka_consensus:
    input:
        basecalls=rules.files.params.reads,
        draft=rules.files.params.ref,
        sam= rules.minimap2.output.sam
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.medaka_consensus.log")
    params:
        outdir=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}"),
        cns_mod = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","cns.mod.fasta")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.sam:q} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft:q} > {params.cns_mod:q}
            medaka_consensus -i {input.basecalls:q} -d {params.cns_mod:q} -o {params.outdir:q} -t 2 &> {log:q}
        else
            touch {output.consensus:q}
        fi
        """

rule gather_cns:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus.fasta"), reference=REFERENCES)
    output:
        os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta")
    run:
        with open(output[0],"w") as fw:
            for in_file in input:
                reference = os.path.basename(os.path.dirname(in_file))
                for record in SeqIO.parse(in_file,"fasta"):
                    fw.write(f">{reference}|{BARCODE}\n{record.seq}\n")