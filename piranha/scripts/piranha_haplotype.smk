

import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.haplo_functions import *
from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *
 

BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],HAPLOTYPING_CONFIG),
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","freebayes.vcf"), reference=REFERENCES)

rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule rasusa:
    input:
        reads= rules.files.params.reads,
        ref=rules.files.params.ref
    params:
        depth = config[KEY_HAPLOTYPE_SAMPLE_SIZE]
    output:
        fastq= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","downsample.fastq")
    run:
        ref_len = 0
        for record in SeqIO.parse(input.ref,"fasta"):
            ref_len = len(record)
        shell("rasusa -i {input.reads:q} -c {params.depth:q} " + f"-g {ref_len}b" + " -o {output.fastq:q}")


rule minimap_for_bam:
    input:
        ref=rules.files.params.ref,
        fastq = rules.rasusa.output.fastq
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2_for_bam.log")
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","downsample.bam")
    shell:
        """
        minimap2 -t {threads} -ax asm20 --secondary=no  \
        {input.ref:q} \
        {input.fastq:q} | samtools sort -@ {threads} -O bam -o {output:q} &> {log:q}
        """

rule freebayes:
    input:
        bam = rules.minimap_for_bam.output.bam,
        ref = rules.files.params.ref
    params:
        allele_freq = config[KEY_MIN_ALLELE_FREQUENCY]
    output:
        vcf = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","freebayes.vcf")
    shell:# where i is for no indels, and X is for no mnp, and u is for complex obs
        """
        freebayes -f {input.ref:q} \
                -F {params.allele_freq} \
                --pooled-continuous \
                -i -X -u \
                {input.bam:q} > {output.vcf:q}
        """

rule flopp:
    input:
        bam = rules.minimap_for_bam.output.bam,
        vcf = rules.freebayes.output.vcf
    threads: workflow.cores
    params:
        ploidy = config[KEY_MAX_HAPLOTYPES],
        partition_path = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_partitions")
    output:
        flopp = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_out.txt"),
        partition = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","haplotyping","flopp_partitions","{reference}_part.txt")
    shell:
        """
        VARIANTS=$(grep -v '#' {input.vcf} | wc -l)
        if [[ $VARIANTS -gt 0 ]]
        then
            flopp -b {input.bam:q} \
            -c {input.vcf:q} \
            -p {params.ploidy} -m \
            -t {threads} \
            -o {output.flopp:q} \
            -P {params.partition_path} 
        else
            echo "No variants called - single reference haplotype"
            touch {output.partition:q}
            touch {output.flopp:q}
        fi
        """

rule haplotype_qc:
    input:
        partition = rules.flopp.output.partition,
        reads = rules.files.params.reads,
        ref = rules.files.params.ref
    params:
        reference = "{reference}",
        min_distance = config[KEY_MIN_HAPLOTYPE_DISTANCE],
        min_reads = config[KEY_MIN_HAPLOTYPE_DEPTH],
        haplodir = os.path.join(config[KEY_TEMPDIR],"reference_groups")
    output:
        txt=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}.haplotypes.txt")
    run:
        ## qc steps
        """
        merge close haplotypes
        evenness statistic, standard deviation
        how many reads includued
        minimum reads for the haplotype to be used
        """
        partitions = parse_partition_file(input.partition)
        seq_index = SeqIO.index(input.reads, "fastq")
        merge_info = collapseClose(rules.flopp.output.flopp,config[KEY_MIN_HAPLOTYPE_DISTANCE],rules.freebayes.output.vcf)

        with open(output.txt,"w") as fhaplo:
            merged_haplo_count = 0
            for subset in merge_info:
                reads = []
                for part in subset:
                    reads.append(partitions[part])
                print(green(f"Haplotype {merged_haplo_count}:"), len(reads))

                if len(reads) > params.min_reads:
                    haplotype = f"{params.reference}.HAP0{merged_haplo_count}"
                    with open(os.path.join(params.haplodir,f"{haplotype}.fastq"),"w") as fw:
                        records = []
                        for read in reads:
                            record = seq_index[read]
                            records.append(record)

                        SeqIO.write(records,fw,"fastq")

                    haplo_ref = os.path.join(params.haplodir,f"{haplotype}.reference.fasta")
                    shell(f"cp {input.ref} {haplo_ref}")
        
                    fhaplo.write(f"{haplotype}\n")
                merged_haplo_count += 1


rule write_yaml:
    input:
        expand(rules.haplotype_qc.output.txt, reference=REFERENCES)
    output:
        yaml = os.path.join(config[KEY_TEMPDIR],HAPLOTYPING_CONFIG)
    run:
        haplotypes = []
        for haplo_file in input:
            with open(haplo_file, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    haplotypes.append(l)
        barcode_config = config
        barcode_config[BARCODE] = haplotypes
        with open(output.yaml, 'w') as fw:
            yaml.dump(barcode_config, fw) 

        print(green(f"Haplotypes for {BARCODE}:"))
        for h in haplotypes:
            print(f"- {h}")
        print(yellow("-----------------------"))

