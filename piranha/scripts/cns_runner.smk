import os
from piranha.utils.config import *
from piranha.analysis.clean_racon import clean_racon_gaps
barcode = config["barcode"]

rule all:
    input:
        expand(os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus.fasta"), taxid=config["haplotypes"])
        # expand(os.path.join(config[KEY_OUTDIR],"consensus_sequences","{taxid}.fasta"), taxid=config["taxa_present"])

rule files:
    params:
        ref=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fasta"),
        reads=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fastq")

rule minimap2_racon0:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref,
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.ref.paf")
    shell:
        """
        minimap2 -x map-ont {input.ref} {input.reads} > {output}
        """


rule racon1:
    input:
        reads=rules.files.params.reads,
        fasta=rules.files.params.ref,
        paf= os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.ref.paf")
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon1.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule clean1:
    input:
        fasta = rules.racon1.output,
        ref = rules.files.params.ref
    params:
        taxon = "{taxid}",
        to_align = os.path.join(config[KEY_TEMPDIR],"racon1_ref.fasta"),
        aligned = os.path.join(config[KEY_TEMPDIR],"racon1_ref.aln.fasta")
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon1.clean.fasta")
    run:
        shell("""
        cat {input.ref} {input.fasta} > {params.to_align:q} && \
        mafft {params.to_align} > {params.aligned}
              """)

        clean_racon_gaps("1", params.taxon, params.aligned, output[0])

rule minimap2_racon1:
    input:
        reads=rules.files.params.reads,
        ref=rules.clean1.output,
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.racon1.paf")
    shell:
        "minimap2 -n 1 -m 1 -M 1 {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads=rules.files.params.reads,
        fasta=rules.clean1.output,
        paf= rules.minimap2_racon1.output
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon2.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule clean2:
    input:
        fasta = rules.racon2.output,
        ref = rules.files.params.ref
    params:
        taxon = "{taxid}",
        to_align = os.path.join(config[KEY_TEMPDIR],"racon2_ref.fasta"),
        aligned = os.path.join(config[KEY_TEMPDIR],"racon2_ref.aln.fasta")
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon2.clean.fasta")
    run:
        shell("""
        cat {input.ref} {input.fasta} > {params.to_align:q} && \
        mafft {params.to_align} > {params.aligned}
              """)

        clean_racon_gaps("2", params.taxon, params.aligned, output[0])

rule minimap2_racon2:
    input:
        reads=rules.files.params.reads,
        ref=rules.clean2.output
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.racon2.paf")
    shell:
        "minimap2 -n 1 -m 1 -M 1 {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads=rules.files.params.reads,
        fasta=rules.clean2.output,
        paf= rules.minimap2_racon2.output
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon3.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule clean3:
    input:
        fasta = rules.racon2.output,
        ref = rules.files.params.ref
    params:
        taxon = "{taxid}",
        to_align = os.path.join(config[KEY_TEMPDIR],"racon3_ref.fasta"),
        aligned = os.path.join(config[KEY_TEMPDIR],"racon3_ref.aln.fasta")
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon3.clean.fasta")
    run:
        shell("""
        cat {input.ref} {input.fasta} > {params.to_align:q} && \
        mafft {params.to_align} > {params.aligned}
              """)

        clean_racon_gaps("3", params.taxon, params.aligned, output[0])


rule minimap2_racon3:
    input:
        reads=rules.files.params.reads,
        ref=rules.clean3.output,
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.racon3.paf")
    shell:
        "minimap2 -n 1 -m 1 {input.ref} {input.reads} > {output}"


rule racon4:
    input:
        reads=rules.files.params.reads,
        fasta=rules.clean3.output,
        paf= rules.minimap2_racon3.output
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon4.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule clean4:
    input:
        fasta = rules.racon3.output,
        ref = rules.files.params.ref
    params:
        taxon = "{taxid}",
        to_align = os.path.join(config[KEY_TEMPDIR],"racon4_ref.fasta"),
        aligned = os.path.join(config[KEY_TEMPDIR],"racon4_ref.aln.fasta")
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon4.clean.fasta")
    run:
        shell("""
        cat {input.ref} {input.fasta} > {params.to_align:q} && \
        mafft {params.to_align} > {params.aligned}
              """)

        clean_racon_gaps("4", params.taxon, params.aligned, output[0])


rule minimap2_racon4:
    input:
        reads=rules.files.params.reads,
        ref=rules.clean4.output,
    output:
        os.path.join(config[KEY_OUTDIR],"{taxid}","racon","mapped.racon4.paf")
    shell:
        "minimap2 -n 1 -m 1 {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        basecalls=rules.files.params.reads,
        draft=rules.clean4.output,
        paf= rules.minimap2_racon4.output
    params:
        outdir=os.path.join(config[KEY_OUTDIR],"{taxid}","medaka"),
        cns_mod = os.path.join(config[KEY_OUTDIR],"{taxid}","racon","racon4.mod.fasta")
    output:
        consensus= os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.paf} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft} > {params.cns_mod}
            medaka_consensus -i {input.basecalls} -d {params.cns_mod} -o {params.outdir} -t 2
        else
            touch {output.consensus}
        fi
        """

# rule gather:
#     input:
#         os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus.fasta")
#     output:
#         os.path.join(config[KEY_OUTDIR],"consensus_sequences","{taxid}.fasta")
#     shell:
#         """
#         cp {input[0]} {output[0]} 
#         """

# rule align_to_reference:
#     input:
#         asm = os.path.join(config[KEY_OUTDIR],"{taxid}","medaka","consensus.fasta"),
#         reference=rules.files.params.ref
#     log: os.path.join(config[KEY_TEMPDIR], "logs/minimap2_sam.log")
#     output:
#         aln = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","consensus.aln.fasta")
#     shell:
#         """
#         minimap2 -a -x asm10 --score-N=0 --sam-hit-only --secondary=no -t  {workflow.cores} \
#         {input.reference:q} '{input.asm}' | \
#         gofasta sam toMultiAlign \
#             -t {workflow.cores} \
#             --reference {input.reference:q} \
#             > '{output.aln}'\
#         """

# rule call_snps:
#     input:
#         reference=rules.files.params.ref,
#         aln = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","consensus.aln.fasta")
#     output:
#         snps = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","snps.csv")
#     shell:
#         """
#         gofasta snps -r -q -o
#         """

# rule identify_excess_snps:
#     input:
#     output:
#     shell:
#         """
#         """

# rule make_snipit_alignment:
#     input:
#         reference=rules.files.params.ref,
#         aln = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","consensus.aln.fasta")
#     output:
#         aln = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","snipit.aln.fasta")
#     run:




# rule make_snipit:
#     input:
#         aln = os.path.join(config[KEY_OUTDIR],"{taxid}","alignment","consensus.aln.fasta")
#     output:
#         svg = os.path.join(config[KEY_OUTDIR],"{taxid}","figures","snipit_plot.svg")
#     shell:
#         """
#         snipit {input.aln:q} -r "outgroup" -o {params.out_stem} -f svg
#         """