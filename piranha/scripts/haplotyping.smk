import os
from piranha.analysis.lengthFilter import gather_filter_reads_by_length 
from piranha.analysis.haplotyping_functions import *
runName = "Pak_ES_Run18"
readSource = './test_data/Pakistan_ES/Run18/'
referenceFile = "references/references.vp1.fasta"
barcodes = []
for filename in os.listdir(readSource):
    barcodes.append(filename)

# first few rules are light read filtering and mapping/getting out what references are found
#probably redundant once integrated into piranha

rule all:
  input: 
    expand(os.path.join(runName,"{barcode}.txt"),barcode=barcodes)

rule filterLength:
    input: 
    params:
        readDir = readSource + "{barcode}"
    output: 
        filtered = readSource + "{barcode}/filtered_reads.fastq"
    run:
        #slight mods from the function in piranha 
        gather_filter_reads_by_length(params.readDir,output.filtered)


rule mapReadsAndGatherRefs:
    input:
        reads = rules.filterLength.output.filtered,
        ref = referenceFile
    params:
        barcode = "{barcode}",
        minReads = 300,
        runName = runName
    output:
        mapped_reads = os.path.join("mappedReads",runName,"{barcode}","mapped_reads.bam"),
        index = os.path.join("mappedReads",runName,"{barcode}","mapped_reads.bai"),
        referencesFound = os.path.join("mappedReads",runName,"{barcode}","referencesFound.txt")
    threads: workflow.cores
    log:
        stdout = os.path.join("logs","mappingStep","minimap",runName,"{barcode}.log"),
        stderr = os.path.join("logs","mappingStep","minimap",runName,"{barcode}.err"),
        samtoolsErr = os.path.join("logs","mappingStep","samtools",runName,"{barcode}.err")
    run:
        #minimap (probably not actually necessary once actually incorporated)
        shell("""
        minimap2 -ax map-ont {input.ref} {input.reads} --secondary=no -o mappedReads/{params.runName}/{params.barcode}/mapped_reads.sam -t {threads} 2> {log.stderr} > {log.stdout}
        samtools view -b -F 4 --threads {threads} mappedReads/{params.runName}/{params.barcode}/mapped_reads.sam | samtools sort --threads {threads} -o mappedReads/{params.runName}/{params.barcode}/mapped_reads.bam 2> {log.samtoolsErr}
        samtools index mappedReads/{params.runName}/{params.barcode}/mapped_reads.bam mappedReads/{params.runName}/{params.barcode}/mapped_reads.bai 2>> {log.samtoolsErr}
        samtools idxstats {output.mapped_reads} | awk '{{if($3>={params.minReads}) print $0}}' > {output.referencesFound} 2>> {log.samtoolsErr}
        """)

checkpoint splitBambyRef:
    input:
        bam = rules.mapReadsAndGatherRefs.output.mapped_reads,
        refs = rules.mapReadsAndGatherRefs.output.referencesFound
    output:
        refBamDir = directory(os.path.join("mappedReads",runName,"{barcode}","referenceBams")),
        referenceFastas = directory(os.path.join("references","{barcode}"))
    params:
        refs = referenceFile
    log:
        stdout = os.path.join("logs","referenceSplitting",runName,"{barcode}.log"),
        stderr = os.path.join("logs","referenceSplitting",runName,"{barcode}.err")
    shell:
        """
        mkdir {output.refBamDir}
        mkdir {output.referenceFastas}
        for reference in $(cut -f1 {input.refs})
        do
        header=$(grep $reference {params.refs})
        sequence=$(grep -A1 $reference {params.refs} | tail -1)
        DN=$(echo $header | cut -d ' ' -f 2 | cut -d '=' -f 2)
        if [[ $DN == "NonPolioEV" ]]
        then
            continue
        fi
        echo $header > {output.referenceFastas}/${{reference}}.fasta
        echo $sequence >> {output.referenceFastas}/${{reference}}.fasta
        samtools view -b --threads {threads} {input.bam} ${{reference}} | samtools sort --threads {threads} -o {output.refBamDir}/${{reference}}.bam > {log.stdout} 2> {log.stderr}
        samtools index {output.refBamDir}/${{reference}}.bam {output.refBamDir}/${{reference}}.bai >> {log.stdout} 2>> {log.stderr}
        done
        """

rule downSampleBam:
    input:
        bam = os.path.join("mappedReads",runName,"{barcode}","referenceBams","{refSeq}.bam")
    output:
        ds_bam = os.path.join("mappedReads",runName,"{barcode}","downsampled","{refSeq}_ds.bam")
    params:
        read_limit = 3000 #would be default or user specified, add later
    log:
        stdout = os.path.join("logs","downsampling",runName,"{barcode}_{refSeq}.err"),
        stderr = os.path.join("logs","downsampling",runName,"{barcode}_{refSeq}.log")
    shell:
        #need to know what proportion to downsample to - so need to get number of mapped reads with flagstat and then get what proportion of that the read limit is
        #I have just written an entire script to make my life easier
        #there is probably a less convoluted way of doing this
        "bash createDownsample.sh {input.bam} {params.read_limit} {output.ds_bam} > {log.stdout} 2> {log.stderr}"

rule variantCall:
    input:
        bam = rules.downSampleBam.output.ds_bam,
        ref = "references/{barcode}/{refSeq}.fasta"
    output:
        vcf = os.path.join("vcfFiles",runName,"{barcode}","{refSeq}.vcf"),
        int_vcf = os.path.join("tmp","vcfFiles",runName,"{barcode}","{refSeq}_int.vcf")
        #need to add in check for zero variants in next rule or flopp will panic
    params:
        alleleFreq = 0.07 #another thing that user should be able to specify
    log:
        stderr = os.path.join("logs","variantCalling",runName,"{barcode}_{refSeq}.err")
    run:
        shell("freebayes -f {input.ref} -F {params.alleleFreq} --pooled-continuous --no-indels --no-mnps --no-complex {input.bam} > {output.int_vcf} 2> {log.stderr}")
        write_contig_headers_vcf(output.int_vcf,output.vcf)


rule callHaplos:
    input:
        bam = rules.downSampleBam.output.ds_bam,
        vcf = os.path.join("vcfFiles",runName,"{barcode}","{refSeq}.vcf")
    params:
        ploidy = 4, #again, needs the option for user to specify
        partitionPath = os.path.join("floppFiles",runName,"{barcode}","partitions")
    output:
        flopp = os.path.join("floppFiles",runName,"{barcode}","{refSeq}.results.flopp"),
        partFile = os.path.join("floppFiles",runName,"{barcode}","partitions","{refSeq}_part.txt")

    log:
        stdout = os.path.join("logs","haplotyping",runName,"{barcode}_{refSeq}.log"),
        stderr = os.path.join("logs","haplotyping",runName,"{barcode}_{refSeq}.err")
    shell:
        """
        VARIANTS=$(grep -v '#' {input.vcf} | wc -l)
        if [[ $VARIANTS -gt 0 ]]
        then
            ~/flopp/target/release/flopp -b {input.bam} -c {input.vcf} -p {params.ploidy} -m -o {output.flopp} -t 1 -P {params.partitionPath} > {log.stdout} 2> {log.stderr}
        else
            echo "No variants called - single reference haplotype" > {output.flopp}
        fi
        """

checkpoint createReadPartition:
    input:
        partFile = rules.callHaplos.output.partFile,
        bam = rules.downSampleBam.output.ds_bam
    output:
        partBamDir = directory(os.path.join("partitionedBams",runName,"{barcode}","{refSeq}"))
    params:
        outPrefix = os.path.join("partitionedBams",runName,"{barcode}","{refSeq}","haplotype")
    run:
        #want to put this script into a function and just import
        #also mod to ignore unmapped
            shell("mkdir {output.partBamDir}")
            get_bam_partition(input.partFile,input.bam,params.outPrefix)
            shell("rm {params.outPrefix}-not_mapped*")

rule createFastqs:
    input:
        partBam = os.path.join("partitionedBams",runName,"{barcode}","{refSeq}","haplotype{haplotype}.bam")
    output:
        partFastq = os.path.join("partitionedReads",runName,"{barcode}","{refSeq}","haplotype{haplotype}.fastq")
    log:
        stderr = os.path.join("logs","samtools_fastq",runName,"{barcode}_{refSeq}_{haplotype}.log")
    shell:
        "samtools fastq {input.partBam} > {output.partFastq} 2> {log.stderr}"


rule runMedaka:
    input:
        fastq = rules.createFastqs.output.partFastq,
        reference = "references/{barcode}/{refSeq}.fasta"
    output:
        medakaVCF = os.path.join("medakaFiles",runName,"{barcode}","{refSeq}","{haplotype}","medaka.vcf"),
    params:
        outDir = os.path.join("medakaFiles",runName,"{barcode}","{refSeq}","{haplotype}")
    log:
        #currently giving an error about files not having the same wildcards?
        stdout = os.path.join("logs","medaka",runName,"{barcode}_{refSeq}_{haplotype}.log"),
        stderr = os.path.join("logs","medaka",runName,"{barcode}_{refSeq}_{haplotype}.err")
    shell:
        """
            medaka_haploid_variant -i {input.fastq} -r {input.reference} -o {params.outDir} -f -x > {log.stdout} 2> {log.stderr}
        """


def aggregate_input_for_alignment(wildcards):
    #Combining all haplotypes into one output per reference seq
    checkpoint_output = checkpoints.createReadPartition.get(**wildcards).output[0]
    return expand(os.path.join("medakaFiles",runName,"{barcode}","{refSeq}","{haplotype}","medaka.vcf"),
           barcode=wildcards.barcode,
           refSeq=wildcards.refSeq,
           haplotype=glob_wildcards(os.path.join(checkpoint_output,"haplotype{haplotype}.bam")).haplotype)

rule createAlignment:
    input:
        #consensus seqs
        variantCalls = aggregate_input_for_alignment,
        ref = "references/{barcode}/{refSeq}.fasta"
    output:
        alignedSeqs = os.path.join("snipitFiles",runName,"{barcode}","{refSeq}","alignedConsensus.fasta")
    run:
        #generates haplo seqs from ref seq and medaka-called variants
        parsedVCF = parseVCFCalls(input.variantCalls)
        makeHaploFasta(parsedVCF,output[0],wildcards.refSeq,input.ref)

            

rule generateSnipit:
    input:
        alignment = rules.createAlignment.output.alignedSeqs
    params:
        ref = "{refSeq}",
        outDir = os.path.join("snipits",runName,"{barcode}","{refSeq}")
    output:
        svg = os.path.join("snipits",runName,"{barcode}","{refSeq}","snp_plot.svg")
    shell:
        "snipit {input.alignment} -r {params.ref} -f svg -d {params.outDir}"

def aggregate_input_for_report(wildcards):
    #combining into one output per barcode
    checkpoint_output = checkpoints.splitBambyRef.get(**wildcards).output[0]
    return (expand(os.path.join("snipits",runName,"{barcode}","{refSeq}","snp_plot.svg"),
           barcode=wildcards.barcode,
           refSeq=glob_wildcards(os.path.join(checkpoint_output,"{refSeq}.bam")).refSeq))

rule generateReport:
    #this is a placeholder rule currently
    input:
        aggregate_input_for_report
    output:
        os.path.join(runName,"{barcode}.txt")
    shell:
        "ls {input} > {output}"