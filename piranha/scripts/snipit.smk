import os
from Bio import SeqIO

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *

rule all:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{taxid}.svg"), taxid=config[KEY_TAXA_PRESENT]),
        os.path.join(config[KEY_TEMPDIR],"snipit.txt")

rule files:
    params:
        ref=os.path.join(config[KEY_TAXA_OUTDIR],"{taxid}.fasta")

rule make_snipit_alignments:
    input:
        refs = config[KEY_REFERENCE_SEQUENCES],
        fasta = config["cns"]
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"snipit","{taxid}.aln.fasta")
    run:
        with open(output.aln,"w") as fw:
            for record in SeqIO.parse(input.ref,"fasta"):
                fw.write(f">{record.id}|reference\n{record.seq}\n")

            for record in SeqIO.parse(input.fasta,"fasta"):
                fw.write(f">{record.id}\n{record.seq}\n")

rule run_snipit:
    input:
        aln = rules.make_snipit_alignments.output.aln
    params:
        out_stem =os.path.join(config[KEY_TEMPDIR],"snipit","{taxid}")
    output:
        os.path.join(config[KEY_TEMPDIR],"snipit","{taxid}.svg")
    run:
        try:
            shell("""snipit {input.aln:q} -o {params.out_stem} -f svg -c wes""")
        except:
            shell("touch {output[0]:q}")

rule gather_graphs:
    input:
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{taxid}.svg"), taxid=config[KEY_TAXA_PRESENT])
    output:
        os.path.join(config[KEY_TEMPDIR],"snipit.txt")
    shell:
        "touch {output:q}"