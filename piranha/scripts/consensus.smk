

import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta"),
        os.path.join(config[KEY_TEMPDIR],"variants.csv"),
        os.path.join(config[KEY_TEMPDIR],"masked_variants.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg"), reference=REFERENCES),

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
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant","consensus.fasta")
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
        """

rule medaka_haploid_variant_cns:
    input:
        reads=rules.files.params.reads,
        ref=rules.medaka_haploid_variant.output.cns
    params:
        model = config[KEY_MEDAKA_MODEL],
        outdir = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns")
    output:
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus_probs.hdf"),
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","medaka.vcf"),
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_haploid_variant_cns","consensus.fasta")
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}_cns.hapoid_variant.log")
    shell:
        """
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
        """


rule join_cns_ref:
    input:
        ref=rules.files.params.ref,
        medaka_cns=rules.medaka_haploid_variant.output.cns,
        cns_cns=rules.medaka_haploid_variant_cns.output.cns
    params:
        reference = "{reference}"
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.fasta")
    run:
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.ref,KEY_FASTA):
                display_name = ""
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]

                fw.write(f">{display_name} {record.description}\n{record.seq}\n")
            if "Sabin" in params.reference:
                for record in SeqIO.parse(input.medaka_cns,KEY_FASTA):
                    record_name = str(SAMPLE).replace(" ","_")
                    fw.write(f">{record_name}\n{record.seq}\n")
            else:
                for record in SeqIO.parse(input.cns_cns,KEY_FASTA):
                    record_name = str(SAMPLE).replace(" ","_")
                    fw.write(f">{record_name}\n{record.seq}\n")


rule align_cns_ref:
    input:
        rules.join_cns_ref.output.fasta
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.aln.fasta")
    shell:
        """
        mafft {input:q} > {output:q}
        """

rule curate_variants:
    input:
        aln = rules.align_cns_ref.output.aln
    params:
        reference = "{reference}"
    output:
        masked = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","masked.csv"),
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns_clean.fasta")
    run:
        masked = clean_medaka_cns(SAMPLE, input.aln, output.fasta)
        with open(output.masked,"w") as fw:
            for var in masked:
                site = int(var) + 1
                fw.write(f"{params.reference},{site},{masked[var]}\n")

rule gather_masked_variants:
    input:
        expand(rules.curate_variants.output.masked, reference=REFERENCES)
    output:
        masked = os.path.join(config[KEY_TEMPDIR],"masked_variants.csv")
    run:
        with open(output.masked,"w") as fw:
            fw.write("reference,site,variant\n")
            for i in input:
                with open(i, "r") as f:
                    for l in f:
                        l=l.rstrip("\n")
                        fw.write(l + '\n')


rule join_clean_cns_ref:
    input:
        ref=rules.files.params.ref,
        cns=rules.curate_variants.output.fasta
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.fasta")
    run:
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.ref,KEY_FASTA):
                display_name = ""
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]

                fw.write(f">{display_name} {record.description}\n{record.seq}\n")
            for record in SeqIO.parse(input.cns,KEY_FASTA):
                record_name = SAMPLE.replace(" ","_")
                fw.write(f">{record_name}\n{record.seq}\n")

rule align_clean_cns_ref:
    input:
        rules.join_clean_cns_ref.output.fasta
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.aln.fasta")
    shell:
        """
        mafft {input:q} > {output:q}
        """

rule make_snipit_graph:
    input:
        aln = rules.align_clean_cns_ref.output.aln
    params:
        out_stem = os.path.join(config[KEY_TEMPDIR],"snipit","{reference}")
    output:
        os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg")
    run:
        try:
            shell("""snipit {input.aln:q} -o {params.out_stem} -f svg -c wes""")
        except:
            shell("touch {output[0]:q}")

rule assess_variants:
    input:
        rules.align_clean_cns_ref.output.aln
    params:
        reference = "{reference}"
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","variants.csv")
    run:
        parse_variants(input[0],output.csv,BARCODE,params.reference)

rule gather_variants:
    input:
        expand(rules.assess_variants.output.csv, reference=REFERENCES)
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"variants.csv")
    run:
        join_variant_files(VARIANT_CALLS_HEADER_FIELDS,input,output.csv)

rule gather_cns:
    input:
        variants = rules.gather_variants.output.csv,
        seqs = expand(rules.curate_variants.output.fasta, reference=REFERENCES)
    output:
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta")
    run:
        var_dict = {}
        with open(input.variants, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                var_dict[row[KEY_REFERENCE]] = row

        with open(output[0],"w") as fw:
            for in_file in input.seqs:
                reference = os.path.basename(os.path.dirname(in_file))
                for record in SeqIO.parse(in_file,KEY_FASTA):
                    variant_count = var_dict[reference]["variant_count"]
                    variant_string = var_dict[reference]["variants"]
                    fw.write(f">{reference}|{BARCODE}|{variant_count}|{variant_string}\n{record.seq}\n")