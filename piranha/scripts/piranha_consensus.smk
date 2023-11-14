#!/usr/bin/env python3


import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"consensus_config.yaml")

rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule medaka_haploid_variant:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    params:
        model = config[KEY_MEDAKA_MODEL],
        outdir =  os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus_probs.hdf"),
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","medaka.vcf"),
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus.fasta"),
        pub_vcf = os.path.join(config[KEY_TEMPDIR],"variant_calls","{reference}.vcf")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.hapoid_variant.log")
    shell:
        """
        [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}
        if [ -s {input.ref:q} ]
        then
            medaka_haploid_variant -i {input.reads:q} \
                                -r {input.ref:q} \
                                -o {params.outdir:q} \
                                -f -x && \
            medaka stitch {output.probs:q} {input.ref:q} {output.cns:q}
        else
            touch {output.cns:q}
            touch {output.probs:q}
            touch {output.vcf:q}
        fi
        cp {output.vcf:q} {output.pub_vcf:q}
        """

rule medaka_haploid_variant_cns:
    input:
        reads=rules.files.params.reads,
        ref=rules.medaka_haploid_variant.output.cns
    params:
        reference = "{reference}",
        model = config[KEY_MEDAKA_MODEL],
        outdir = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus_probs.hdf"),
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","medaka.vcf"),
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus.fasta")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}_cns.hapoid_variant.log")
    run:
        if not "Sabin" in params.reference:
            shell("""
            [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}
            if [ -s {input.ref:q} ]
            then
                medaka_haploid_variant -i {input.reads:q} \
                                    -r {input.ref:q} \
                                    -o {params.outdir} \
                                    -f -x && \
                medaka stitch {output.probs} {input.ref} {output.cns}
            else
                touch {output.cns:q}
                touch {output.probs:q}
                touch {output.vcf:q}
            fi
            """)
        else:
            shell("""
            touch {output.cns:q}
            touch {output.probs:q}
            touch {output.vcf:q}
            """)

rule gather_merge_cns:
    input:
        ref=expand(rules.files.params.ref,reference=REFERENCES),
        cns = expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus.fasta"),reference=REFERENCES)
        # cns_cns = expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus.fasta"),reference=REFERENCES)
    output:
        yaml = os.path.join(config[KEY_TEMPDIR],"consensus_config.yaml")
    run:
        sequences = collections.defaultdict(set)
        ref_seqs = {}
        for cns_file in input.cns:
            haplodir = "/".join(cns_file.split("/")[:-1])
            haplo_bam = os.path.join(haplodir, "calls_to_ref.bam")
            for record in SeqIO.parse(cns_file, "fasta"):
                print(record.id)
                sequences[str(record.seq)].add(haplo_bam)
                ref_seqs[str(record.seq)] = record.id
        
        ref_file_dict = {}
        for ref_file in input.ref:
            for record in SeqIO.parse(ref_file,"fasta"):
                ref_file_dict[record.id] = ref_file

        for seq in sequences:
            print(seq[:5],"...",sequences[seq])
        
        ref_count = collections.Counter()
        ref_cns = []
        new_config = config

        for seq in sequences:
            ref = ref_seqs[seq]
            ref_count[ref] +=1
            record_id = f"{ref_seqs[seq]}.CNS00{ref_count[ref]}"
            ref_cns.append(record_id)

            haplo_bams = ""
            for i in sequences[seq]:
                haplo_bams += f"'{i}' "

            out_bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{record_id}.merged_cns.bam")
            read_count_out = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{record_id}.read_count.txt")
            shell(f"samtools merge -f '{out_bam}' {haplo_bams} && samtools index '{out_bam}'")
            shell(f"samtools view -F 0x904 -c '{out_bam}' > {read_count_out}")

            with open(read_count_out,"r") as f:
                for l in f:
                    read_count = l.rstrip()
                    new_config[record_id] = read_count

            with open(os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{record_id}.merged_cns.fasta"),"w") as fseq:
                fseq.write(f">{record_id} read_count={read_count}\n{seq}\n")

            cns_ref = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{record_id}.ref.fasta")
            shell(f"cp '{ref_file_dict[ref]}' '{cns_ref}'")
            
            
        new_config[BARCODE] = ref_cns

        with open(output.yaml, 'w') as fw:
            yaml.dump(new_config, fw) 
        print(yellow("-----------------------"))


