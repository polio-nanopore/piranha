import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

BARCODE = config[KEY_BARCODE]
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta"),
        os.path.join(config[KEY_OUTDIR],"variants.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","variants.vcf"), reference=REFERENCES),
        os.path.join(config[KEY_OUTDIR],"variation_info.json"),
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg"), reference=REFERENCES)

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
    output:
        cns_mod = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","cns.mod.fasta"),
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.sam:q} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft:q} > {output.cns_mod:q}
            medaka_consensus -i {input.basecalls:q} -d {output.cns_mod:q} -o {params.outdir:q} -t 2 &> {log:q}
        else
            touch {output.consensus:q}
        fi
        """

rule medaka_variant:
    input:
        ref = rules.medaka_consensus.output.cns_mod,
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus_probs.hdf"),
        bam = rules.sort_index.output.bam
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.medaka_variant.log")
    params:
        temp_vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","variants.pre.vcf")
    output:
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","variants.vcf")
    shell:
        """
        medaka variant {input.ref:q} {input.probs:q} {params.temp_vcf:q}  &> {log} 
        medaka tools annotate {params.temp_vcf:q} {input.ref:q} {input.bam:q} {output[0]} &> {log}
        """

rule sam_to_seq:
    input:
        sam = rules.minimap2.output.sam,
        ref = rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","pseudoaln.fasta")
    shell:
        """
        gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q} &> {log}
        """

rule get_variation_info:
    input:
        expand(rules.files.params.reads, reference=REFERENCES),
        expand(rules.sam_to_seq.output.fasta,reference=REFERENCES)
    output:
        json = os.path.join(config[KEY_OUTDIR],"variation_info.json")
    run:
        # this is for making a figure
        variation_dict = {}
        for reference in REFERENCES:
            ref = os.path.join(config[KEY_TEMPDIR],"reference_groups",f"{reference}.reference.fasta")
            fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","pseudoaln.fasta")
            
            var_dict = get_variation_pcent(ref,fasta)
            variation_dict[reference] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))

rule join_cns_ref:
    input:
        ref=rules.files.params.ref,
        cns=rules.medaka_consensus.output.consensus
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.fasta")
    run:
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.ref,"fasta"):
                fw.write(f">{record.description}\n{record.seq}\n")
            for record in SeqIO.parse(input.cns,"fasta"):
                fw.write(f">{BARCODE}\n{record.seq}\n")

rule align_cns_ref:
    input:
        rules.join_cns_ref.output.fasta
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.aln.fasta")
    shell:
        """
        mafft {input:q} > {output:q}
        """

rule make_snipit_graph:
    input:
        aln = rules.align_cns_ref.output.aln
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
        rules.align_cns_ref.output.aln
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
        csv = os.path.join(config[KEY_OUTDIR],"variants.csv")
    run:
        join_variant_files(VARIANT_CALLS_HEADER_FIELDS,input,output.csv)

rule gather_cns:
    input:
        variants = os.path.join(config[KEY_OUTDIR],"variants.csv"),
        seqs = expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","consensus.fasta"), reference=REFERENCES)
    output:
        os.path.join(config[KEY_OUTDIR],"consensus_sequences.fasta")
    run:
        var_dict = {}
        with open(input.variants, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                var_dict[row[KEY_REFERENCE]] = row

        with open(output[0],"w") as fw:
            for in_file in input.seqs:
                reference = os.path.basename(os.path.dirname(in_file))
                for record in SeqIO.parse(in_file,"fasta"):
                    variant_count = var_dict[reference]["variant_count"]
                    variant_string = var_dict[reference]["variants"]
                    fw.write(f">{reference}|{BARCODE}|{variant_count}|{variant_string}\n{record.seq}\n")