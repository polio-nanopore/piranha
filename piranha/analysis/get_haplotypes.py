import os
import csv
from Bio import SeqIO
import collections

from piranha.utils.log_colours import green

def get_variation_pcent(ref,fasta):
    ref_record = SeqIO.read(ref,"fasta")
    
    #set up site counter, 0-based
    variant_sites = {}
    for i in range(len(ref_record)):
        variant_sites[i] = 0

    c = 0
    for record in SeqIO.parse(fasta,"fasta"):
        c +=1
        index = 0
        for pos in zip(str(record.seq),str(ref_record.seq)):
            
            if pos[0] != pos[1]:
                variant_sites[index]+=1
            index +=1
    x = []
    y = []
    for site in variant_sites:
        pcent_variants = round(100*(variant_sites[site]/c), 1)
        x.append(site+1)
        y.append(pcent_variants)
    var_dict = {"x":x,"y":y}
    return var_dict

def parse_vcf(fasta,vcf,min_reads,haplotypes_out):
    variant_sites = []
    with open(vcf,"r") as f:
        for l in f:
            if not l.startswith("#"):
                l = l.rstrip("\n")
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = l.split("\t")
                variant_sites.append(int(POS))
                
    site_str = ';'.join([str(i) for i in variant_sites])
    print(site_str)
    
    haplo_counter = collections.defaultdict(list)
    for record in SeqIO.parse(fasta,"fasta"):
        read_haplotype = "".join([record.seq[i] for i in variant_sites])
        haplo_counter[read_haplotype].append(record.id)
    
    read_haplotypes = {}
    print(green("Haplotypes identified\n--------"))
    with open(haplotypes_out,"w") as fw:
        fw.write(f"sites,haplotype,num_reads,make_cns,read_ids\n")
        for i in haplo_counter:
            num_reads = len(haplo_counter[i])
            read_str = ';'.join(haplo_counter[i])
            if num_reads > min_reads:
                print(green(f"{i}: ") + f"{len(haplo_counter[i])}")
                read_haplotypes[i] = haplo_counter[i]
                fw.write(f"{site_str},{i},{num_reads},True,{read_str}\n")
            else:
                fw.write(f"{site_str},{i},{num_reads},False,{read_str}\n")
    return read_haplotypes

                
def write_haplotype_fastq(reads,read_haplotypes,outdir):
    
    reads = SeqIO.index(reads,"fastq")
    for h in read_haplotypes:
        with open(os.path.join(outdir, f"{h}.fastq"),"w") as fw:
            for read in read_haplotypes[read]:
                record = reads[read]
                SeqIO.write(fw,record,"fastq")
            
def get_haplotypes(fasta,vcf,reads,out_haplotypes,outdir,min_reads):

    haps = parse_vcf(fasta,vcf,min_reads,out_haplotypes)

    write_haplotype_fastq(reads,haps,outdir)

    return list(haps.keys())
