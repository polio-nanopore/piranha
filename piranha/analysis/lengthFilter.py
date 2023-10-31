from Bio import SeqIO
import os
import gzip
import sys

#taken from piranha.analysis.preprocessing
def gather_filter_reads_by_length(dir_in,reads_out):

    # raise Exception(f"I'm a test exception for {dir_in}")

    if not os.path.exists(dir_in):
        os.mkdir(dir_in)

    fastq_records = []

    total_reads = 0
    with open(reads_out,"w") as fw:
        #root, dirs, files
        for r,d,f in os.walk(dir_in):
            for reads_in in f:
                if reads_in.endswith(".gz") or reads_in.endswith(".gzip"):
                    with gzip.open(os.path.join(dir_in,reads_in), "rt") as handle:
                        for record in SeqIO.parse(handle, "fastq"):
                            total_reads +=1
                            length = len(record)
                            if length > 1000 and length < 1300:
                                fastq_records.append(record)
                                

                elif reads_in.endswith(".fastq") or reads_in.endswith(".fq"):
                    for record in SeqIO.parse(os.path.join(dir_in,reads_in),"fastq"):
                        total_reads +=1
                        length = len(record)
                        if length > 1000 and length < 1300:
                            fastq_records.append(record)

        # print(green(f"Total reads {barcode}:"),total_reads)
        # print(green(f"Total passed reads {barcode}:"),len(fastq_records))
        SeqIO.write(fastq_records,fw, "fastq")

if __name__ == '__main__':
    #args are barcode dir, output file
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
    
    gather_filter_reads_by_length(snakemake.params.readDir,snakemake.output)