import os
import csv
from Bio import SeqIO
import collections
import sys
import yaml
csv.field_size_limit(sys.maxsize)

from piranha.utils.log_colours import green,cyan

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

    variant_info = {}
    for site in variant_sites:
        variant_info[site] = collections.Counter()

    for record in SeqIO.parse(fasta,"fasta"):
        for site in variant_sites:
            variant = str(record.seq)[site]
            variant_info[site][variant]+=1

    x = []
    y = []
    info = []
    variation_info = []
    for site in variant_sites:
        
        pcent_variants = round(100*(variant_sites[site]/c), 1)
        x = site+1
        y = pcent_variants

        var_counts = variant_info[site]
        site_data = {"Position":x,"Percentage":y}
        for i in var_counts:
            site_data[f"{i} reads"] = var_counts[i]
        
        variation_info.append(site_data)

    return variation_info

def parse_vcf(fasta,vcf,min_reads,min_pcent,taxon,haplotypes_out):
    variant_sites = []
    with open(vcf,"r") as f:
        for l in f:
            if not l.startswith("#"):
                l = l.rstrip("\n")
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = l.split("\t")
                
                variant_sites.append(int(POS))
                
    site_str = ';'.join([str(i) for i in variant_sites])
    
    print(site_str)
    haplo_counter = collections.Counter()
    haplo_records = collections.defaultdict(list)
    read_count = 0
    for record in SeqIO.parse(fasta,"fasta"):
        read_count +=1
        read_haplotype = ""
        for i in variant_sites:
            read_haplotype+=f"{i}{record.seq[i-1]}"
        haplo_counter[read_haplotype] +=1
    print(haplo_counter)
    for record in SeqIO.parse(fasta,"fasta"):
        read_haplotype = ""
        for i in variant_sites:
            read_haplotype+=f"{i}{record.seq[i-1]}"
        if 100*(haplo_counter[read_haplotype]/read_count) > min_pcent:
            haplo_records[read_haplotype].append(record.id)

    read_haplotypes = {}
    print(green("Haplotypes identified\n--------\n--------"))
    print("Sites: ", site_str)
    with open(haplotypes_out,"w") as fw:
        fw.write(f"taxon,sites,haplotype,num_reads,make_cns\n")
        for i in haplo_records:
            num_reads = len(haplo_records[i])
            if num_reads > min_reads:
                print(green(f"{i}: ") + f"{len(haplo_records[i])}")
                read_haplotypes[i] = haplo_records[i]
                fw.write(f"{taxon},{site_str},{i},{num_reads},True\n")
            else:
                fw.write(f"{taxon},{site_str},{i},{num_reads},False\n")
    return read_haplotypes

def write_haplotype_ref(ref,taxon,read_haplotypes,outdir):
    reference = ""
    for record in SeqIO.parse(ref, "fasta"):
        reference = record

    for h in read_haplotypes:
        with open(os.path.join(outdir, f"{taxon}_{h}.reference.fasta"),"w") as fw:
            fw.write(f">{taxon}|reference\n{record.seq}\n")

def write_haplotype_fastq(reads,taxon,read_haplotypes,outdir):
    
    reads = SeqIO.index(reads,"fastq")
    for h in read_haplotypes:
        with open(os.path.join(outdir, f"{taxon}_{h}.fastq"),"w") as fw:
            for read in read_haplotypes[h]:
                record = reads[read]
                SeqIO.write(record,fw,"fastq")


def get_haplotypes(fasta,vcf,reads,ref,out_haplotypes,outdir,taxon,min_reads,min_pcent):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    haps = parse_vcf(fasta,vcf,min_reads,min_pcent,taxon,out_haplotypes)

    write_haplotype_fastq(reads,taxon,haps,outdir)
    write_haplotype_ref(ref,taxon,haps,outdir)

    return list(haps.keys())


def gather_haplotype_data(input_list,output_csv,config_out,config):

    with open(output_csv, "w") as fw:
        writer = csv.DictWriter(fw, fieldnames = ["taxon","sites","haplotype","num_reads","make_cns","read_ids"],lineterminator='\n')
        writer.writeheader()

        for haplotype_file in input_list:

            with open(haplotype_file,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    writer.writerow(row)
    
    with open(config_out,"w") as fw:

        haplotypes_out = []
        with open(output_csv,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row["make_cns"] == "True":
                        taxon = row["taxon"]
                        haplotype = row["haplotype"]
                        haplotype_stem = f"{taxon}_{haplotype}"
                        haplotypes_out.append(haplotype_stem)
        config["haplotypes"] = haplotypes_out
        
        yaml.dump(config, fw)
