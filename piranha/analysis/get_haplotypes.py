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
    x = []
    y = []
    for site in variant_sites:
        pcent_variants = round(100*(variant_sites[site]/c), 1)
        x.append(site+1)
        y.append(pcent_variants)
    var_dict = {"x":x,"y":y}
    return var_dict

def parse_vcf(fasta,vcf,min_reads,min_pcent,taxon,haplotypes_out):
    variant_sites = []
    with open(vcf,"r") as f:
        for l in f:
            if not l.startswith("#"):
                l = l.rstrip("\n")
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = l.split("\t")
                
                variant_sites.append(int(POS))
                
    site_str = ';'.join([str(i) for i in variant_sites])
    
    
    haplo_counter = collections.Counter()
    haplo_records = collections.defaultdict(list)
    read_count = 0
    for record in SeqIO.parse(fasta,"fasta"):
        read_count +=1
        read_haplotype = "".join([record.seq[i] for i in variant_sites])
        haplo_counter[read_haplotype] +=1

    for record in SeqIO.parse(fasta,"fasta"):
        read_haplotype = "".join([record.seq[i] for i in variant_sites])
        if 100*(haplo_counter[read_haplotype]/read_count) > min_pcent:
            haplo_records[read_haplotype].append(record.id)
    
    read_haplotypes = {}
    print(green("Haplotypes identified\n--------\n--------"))
    print("Sites: ", site_str)
    with open(haplotypes_out,"w") as fw:
        fw.write(f"taxon,sites,haplotype,num_reads,make_cns,read_ids\n")
        for i in haplo_records:
            num_reads = len(haplo_records[i])
            read_str = ';'.join(haplo_records[i])
            if num_reads > min_reads:
                print(green(f"{i}: ") + f"{len(haplo_records[i])}")
                read_haplotypes[i] = haplo_records[i]
                fw.write(f"{taxon},{site_str},{i},{num_reads},True,{read_str}\n")
            else:
                fw.write(f"{taxon},{site_str},{i},{num_reads},False,{read_str}\n")
    return read_haplotypes

                
def write_haplotype_fastq(reads,taxon,read_haplotypes,outdir):
    
    reads = SeqIO.index(reads,"fastq")
    for h in read_haplotypes:
        with open(os.path.join(outdir, f"{taxon}_{h}.fastq"),"w") as fw:
            for read in read_haplotypes[h]:
                record = reads[read]
                SeqIO.write(record,fw,"fastq")
            
def get_haplotypes(fasta,vcf,reads,out_haplotypes,outdir,taxon,min_reads,min_pcent):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    haps = parse_vcf(fasta,vcf,min_reads,min_pcent,taxon,out_haplotypes)

    write_haplotype_fastq(reads,taxon,haps,outdir)

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

