import csv
from Bio import SeqIO
import collections

ref = SeqIO.read("analysis_2021-12-14/barcode01/categorised_sample/binned_reads/Sabin2_vacc.fasta","fasta")

variant_sites = {}
for i in range(len(ref)):
    variant_sites[i] = 0
    
c = 0
for record in SeqIO.parse("analysis_2021-12-14/barcode01/categorised_sample/test_pseudoalign.fasta","fasta"):
#     print(record.id)
    c +=1
    for index in range(len(ref)):
        
        if record.seq[index] != ref.seq[index]:
            variant_sites[index]+=1

print(len(ref), len(variant_sites),c)
for site in sorted(variant_sites):
    print(site, variant_sites[site], round(100*(variant_sites[site]/c), 2))

cns = SeqIO.read("analysis_2021-12-14/barcode01/categorised_sample/consensus_sequences/Sabin2_vacc/medaka/consensus.fasta","fasta")
print(cns)

for index in range(len(cns)):
        
    if cns.seq[index] != ref.seq[index]:
        print(index, cns.seq[index],ref.seq[index])
