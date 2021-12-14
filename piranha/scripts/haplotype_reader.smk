import os
import collections
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

rule all:
    input:
        expand(os.path.join(config[KEY_OUTDIR],"{taxid}","haplotype","pseudo_aln.fasta"), taxid=config["taxa_present"])

rule files:
    params:
        ref=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fasta"),
        reads=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fastq")

rule rule_get_read_aln:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    log: os.path.join(config[TEMPDIR],"{taxid}","haplotype.log")
    params:
        sam = os.path.join(config[KEY_OUTDIR],"{taxid}","haplotype","mapped.ref.sam")
    output:
        aln = os.path.join(config[KEY_OUTDIR],"{taxid}","haplotype","pseudo_aln.fasta")
    shell:
        """
        minimap2 -a -x map-ont --sam-hit-only --secondary=no -t {workflow.cores} \
                 {input.ref:q} {input.reads:q} -o {params.sam:q} &> {log:q} 
                    gofasta sam toMultiAlign \
                        -s {params.sam:q} \
                        -t {workflow.cores} \
                        --reference {input.ref:q} \
                         > '{output.aln}'
        """

