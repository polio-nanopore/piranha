#!/usr/bin/env python3
import os
import collections
import csv
from Bio import SeqIO
import yaml

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *
##### Target rules #####

clusters = config[KEY_CLUSTERS]

rule all:
    input:
        expand(os.path.join(config[KEY_OUTDIR],"phylogenetics","{reference_group}.tree"), reference_group=clusters)

rule mafft:
    input:
        seqs = os.path.join(config[KEY_OUTDIR],"phylogenetics","{reference_group}.fasta")
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"phylogenetics","{reference_group}.aln.fasta"),
        pub_aln = os.path.join(config[KEY_OUTDIR],"phylogenetics","{reference_group}.alignment.fasta")
    shell:
        """
        mafft --quiet {input.seqs:q} > {output.aln:q} && cp {output.aln:q} {output.pub_aln:q}
        """

rule iqtree:
    input:
        aln = rules.mafft.output.aln
    output:
        tree = os.path.join(config[KEY_TEMPDIR],"phylogenetics","{reference_group}.aln.fasta.treefile")
    log:
        os.path.join(config[KEY_TEMPDIR],"phylogenetics","{reference_group}.log")
    shell:
        """
        iqtree  -s {input.aln:q} \
                -m HKY \
                -czb \
                -blmin  0.0000000001 \
                -nt 1 \
                -redo \
                --fast \
                -o outgroup \
                -quiet &> {log:q}
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree
    output:
        tree = os.path.join(config[KEY_TEMPDIR],"phylogenetics","{reference_group}.unpruned.tree")
    shell:
        """
        jclusterfunk prune  -i {input.tree:q} \
                            -o {output.tree:q} \
                            -t outgroup \
                            -f newick 
        """

rule annotate:
    input:
        tree = rules.prune_outgroup.output.tree,
        csv = os.path.join(config[KEY_OUTDIR],"phylogenetics", "annotations.csv")
    params:
        annotations = config[KEY_TREE_ANNOTATIONS]
    output:
        tree = os.path.join(config[KEY_OUTDIR],"phylogenetics","{reference_group}.tree")
    shell:
        """
        jclusterfunk annotate -c name \
                             -i {input.tree:q} \
                             -m {input.csv:q} \
                             --tip-attributes {params.annotations}\
                             -f nexus \
                             -o {output.tree:q} \
                            --ignore-missing
        """
