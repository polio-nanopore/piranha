import os

barcode = config["barcode"]

rule all:
    input:
        expand(os.path.join(config["outdir"],"{taxid}","medaka","consensus.fasta"), taxid=config["taxa"]),
        expand(os.path.join(config["outdir"],"consensus_sequences","{taxid}.fasta"), taxid=config["taxa"])

rule minimap2_racon0:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        ref=os.path.join(config["reference_path"],"{taxid}.fasta"),
    output:
        os.path.join(config["outdir"],"{taxid}","racon","mapped.ref.paf")
    shell:
        """
        minimap2 -x map-ont {input.ref} {input.reads} > {output}
        """


rule racon1:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        fasta=os.path.join(config["reference_path"],"{taxid}.fasta"),
        paf= os.path.join(config["outdir"],"{taxid}","racon","mapped.ref.paf")
    output:
        os.path.join(config["outdir"],"{taxid}","racon","racon1.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule minimap2_racon1:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        ref=os.path.join(config["outdir"],"{taxid}","racon","racon1.fasta"),
    output:
        os.path.join(config["outdir"],"{taxid}","racon","mapped.racon1.paf")
    shell:
        "minimap2 -n 1 -m 1 -M 1 {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        fasta=os.path.join(config["outdir"],"{taxid}","racon","racon1.fasta"),
        paf= os.path.join(config["outdir"],"{taxid}","racon","mapped.racon1.paf")
    output:
        os.path.join(config["outdir"],"{taxid}","racon","racon2.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule minimap2_racon2:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        ref=os.path.join(config["outdir"],"{taxid}","racon","racon2.fasta")
    output:
        os.path.join(config["outdir"],"{taxid}","racon","mapped.racon2.paf")
    shell:
        "minimap2 -n 1 -m 1 {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        fasta=os.path.join(config["outdir"],"{taxid}","racon","racon2.fasta"),
        paf= os.path.join(config["outdir"],"{taxid}","racon","mapped.racon2.paf")
    output:
        os.path.join(config["outdir"],"{taxid}","racon","racon3.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule minimap2_racon3:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        ref=os.path.join(config["outdir"],"{taxid}","racon","racon3.fasta"),
    output:
        os.path.join(config["outdir"],"{taxid}","racon","mapped.racon3.paf")
    shell:
        "minimap2 -n 1 -m 1 {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        fasta=os.path.join(config["outdir"],"{taxid}","racon","racon3.fasta"),
        paf= os.path.join(config["outdir"],"{taxid}","racon","mapped.racon3.paf")
    output:
        os.path.join(config["outdir"],"{taxid}","racon","racon4.fasta")
    shell:
        """
        if [ -s {input.paf} ]
        then
            racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}
        else
            touch {output}
        fi
        """

rule minimap2_racon4:
    input:
        reads=os.path.join(config["read_path"],"{taxid}.fastq"),
        ref=os.path.join(config["outdir"],"{taxid}","racon","racon4.fasta"),
    output:
        os.path.join(config["outdir"],"{taxid}","racon","mapped.racon4.paf")
    shell:
        "minimap2 -n 1 -m 1 {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        basecalls=os.path.join(config["read_path"],"{taxid}.fastq"),
        draft=os.path.join(config["outdir"],"{taxid}","racon","racon4.fasta"),
        paf= os.path.join(config["outdir"],"{taxid}","racon","mapped.racon4.paf")
    params:
        outdir=os.path.join(config["outdir"],"{taxid}","medaka"),
        cns_mod = os.path.join(config["outdir"],"{taxid}","racon","racon4.mod.fasta")
    output:
        consensus= os.path.join(config["outdir"],"{taxid}","medaka","consensus.fasta")
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

rule gather:
    input:
        os.path.join(config["outdir"],"{taxid}","medaka","consensus.fasta")
    output:
        os.path.join(config["outdir"],"consensus_sequences","{taxid}.fasta")
    shell:
        """
        cp {input[0]} {output[0]} 
        """
