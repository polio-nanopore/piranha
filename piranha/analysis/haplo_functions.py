#!/usr/bin/env python3
import csv
from Bio import SeqIO
import collections
from piranha.utils.config import *
import os
from scipy.stats import chisquare
import numpy as np

from piranha.utils.log_colours import green,cyan,red


def parse_partition_file(partition_file):

    partitions = collections.defaultdict(set)
    with open(partition_file, "r") as f:
        part = ""
        for l in f:
            l=l.rstrip("\n")
            if l[0] == "#":
                part = l.lstrip("#")
                part = int(part)
            else:
                read_id = l.split("\t")[0]
                partitions[part].add(read_id)

    return partitions



def collapse_close(flopp_file,distance_cutoff,vcf_file):
    haplo_pos = 0
    collapsed = []
    output = []
    final_haplos = []
    vcf = parse_VCF(vcf_file)
    parsed_flopp = parse_flopp_reads(flopp_file)
    for h in parsed_flopp:
        collapse_into_current = [haplo_pos]
        if haplo_pos in collapsed:
            #don't consider haplotypes already merged into another
            haplo_pos += 1
            continue
        new_h = {}
        found_others = False
        
        #for each haplo, check if any others are identical/within distance cutoff, collapse in and combine reads if so, remove collapsed haplo from further consideration
        haplo_pos2 = haplo_pos +1
        for h2 in parsed_flopp[haplo_pos2:]:
            if haplo_pos2 in collapsed:
                #already assigned to a merge group
                continue
            distance = get_haplo_dist(h,h2,vcf)
            if distance <= distance_cutoff and new_h == {}:
                #collapse values and create new haplo 
                for pos in h:
                    #combine
                    new_h[pos] = {}
                    new_h[pos]['total assigned'] = h[pos]['total assigned'] + h2[pos]['total assigned']
                    new_h[pos]['base support counts'] = sum_base_support(h[pos]['base support counts'],h2[pos]['base support counts'])
                    #check if allele assignment is different now that reads are collapsed
                    new_h[pos]['assigned allele'] = max(new_h[pos]['base support counts'],key=lambda x:new_h[pos]['base support counts'][x])
                    
            
                collapsed.append(haplo_pos2)
                collapse_into_current.append(haplo_pos2)
                found_others = True

            elif distance <= distance_cutoff:
                for pos in h:
                    new_h[pos]['total assigned'] = new_h[pos]['total assigned'] + h2[pos]['total assigned']
                    new_h[pos]['base support counts'] = sum_base_support(new_h[pos]['base support counts'],h2[pos]['base support counts'])
                    #check if allele assignment is different now that reads are collapsed
                    new_h[pos]['assigned allele'] = max(new_h[pos]['base support counts'],key=lambda x:new_h[pos]['base support counts'][x])
                        
                collapsed.append(haplo_pos2)
                collapse_into_current.append(haplo_pos2)
                
            haplo_pos2 += 1

        if not found_others:
            new_h = h
        
        #need these merged versions of parsed flopp output to generate even-ness, maybe make one do-all function
        final_haplos.append(new_h)
        haplo_pos += 1
        output.append(collapse_into_current)
        
    return output

def parse_VCF(vcf_file):
    vcf_info = {}
    with open(vcf_file) as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                line = line.strip('\n').split('\t')
                pos, ref, alt, info = line[1], line[3], line[4], line[7]
                pos = int(pos)
                info_dict = dict([x.split('=') for x in info.split(';')])
                vcf_info[pos] = {}
                vcf_info[pos]['0'] = ref
                if ',' in alt:
                    #multiallelic SNP
                    vcf_info[pos]['1'], vcf_info[pos]['2'] = alt.split(',')
                else:
                    vcf_info[pos]['1'] = alt 
                vcf_info[pos]['read_counts'] = {}
                ref_reads = info_dict['RO']
                alt_reads = info_dict['AO']
                #float because freebayes can have 'fractional' observations
                vcf_info[pos]['read_counts']['0'] = float(ref_reads)
                if ',' in alt_reads:
                    counts = [float(x) for x in alt_reads.split(',')]
                    vcf_info[pos]['read_counts']['1'] = counts[0]
                    vcf_info[pos]['read_counts']['2'] = counts[1]
                else:
                    vcf_info[pos]['read_counts']['1'] = float(alt_reads)
                vcf_info[pos]['frequencies'] = {}
                for allele in vcf_info[pos]['read_counts']:
                    vcf_info[pos]['frequencies'][allele] = vcf_info[pos]['read_counts'][allele]/sum(vcf_info[pos]['read_counts'].values())
                #frequencies list allele freqs in 0,1,2 order

    return vcf_info

def parse_flopp_reads(floppPath):
    #returns an array of haplotypes; each haplo is a dict with genomic positions as keys, each pos has a dict with allele, read support
    #might need to make some changes for more than 2 alt alleles
    #TODO change to recording values for all bases
    haps = []
    with open(floppPath) as file:
        file.readline() #header, gives ref mapped to
        for line in file:
            if line.startswith('1:'):
                #first haplo info line, get haplo count and init haplo entries
                line = [x for x in line.strip('\n').split('\t') if x != '']
                pos = line[0].split(':')[1]
                pos=int(pos)
                k = int((len(line) - 1)/2) #number of haplotypes
                hap_cols = line[1:k+1]
                support_cols = line[k+1:]
                
                for i in range(k):
                    haps.append({})
                    out_dict = {}
                    out_dict['assigned allele'] = hap_cols[i]
                    #something has to change here, want reads for all possible alleles
                    out_dict['total assigned'], out_dict['base support counts'] = parse_support_col(support_cols[i])
                    haps[i][pos] = out_dict

            elif line == '*****\n':
                #EOF
                continue

            else:
                line = line.strip('\n').split('\t')
                pos = line[0].split(':')[1]
                pos=int(pos)
                hap_cols = line[1:k+1]
                support_cols = line[k+1:]

                for i in range(k):
                    out_dict = {}
                    out_dict['assigned allele'] = hap_cols[i]
                    out_dict['total assigned'], out_dict['base support counts'] = parse_support_col(support_cols[i])
                    haps[i][pos] = out_dict
    return haps

def parse_support_col(col):
    base_split = col.split('|')
    total_reads = 0
    base_support = {}
    for base in base_split:
        base_val,support_reads = base.split(':')
        # if base_val == selbase_val:
        #     base_support = int(support_reads)
        base_support[base_val] = int(support_reads)
        total_reads += int(support_reads)
    return (total_reads,base_support)

def check_site_sig(base_support,allele_freqs):
    #check if allele freqs are consistent with random sample from overall reads for this site
    #need to know these are the same order
    #reads might have 100% one allele, need to insert zero in this case, then sort to make sure order correct
    for key in allele_freqs:
        if key not in base_support:
            base_support[key] = 0
    read_counts = [base_support[k] for k in sorted(base_support)]
    exp = [allele_freqs[allele_freq]*sum(read_counts) for allele_freq in allele_freqs]
    try:
        pval = chisquare(f_obs=read_counts,f_exp=exp).pvalue
    except Exception as e:
        print('Observed values for read counts',read_counts)
        print('Expected values for read counts',exp)
        print('Allele frequencies',allele_freqs)
        raise ValueError(f'chi squared calculation failed with error, {e}, see above values')
    if pval < 0.01:
        return True
    else:
        return False

def get_haplo_dist(h1,h2,parsed_VCF):
    #returns the number of reliable SNPs separating 2 haplotypes
    sig_pos1, sig_pos2 = [],[]
    # for position in each haplo, check significance (i.e. are the base proportions differnt from a random sample of total read set), add to list if is
    #need a lookup from vcf parsing to get allele_freqs
    for pos in h1:
        if check_site_sig(h1[pos]['base support counts'],parsed_VCF[pos]['frequencies']) and h1[pos]['assigned allele'] != '0':
            #need to append both position and allele to list
            sig_pos1.append((pos,h1[pos]['assigned allele']))
    for pos in h2:
        if check_site_sig(h2[pos]['base support counts'],parsed_VCF[pos]['frequencies']) and h2[pos]['assigned allele'] != '0':
            sig_pos2.append((pos,h2[pos]['assigned allele']))
    #how many of the reliable SNPs are only in one
    #set arith?
    h1_only = set(sig_pos1) - set(sig_pos2)
    h2_only = set(sig_pos2) - set(sig_pos1)
    #this counts the same site with different alt alleles as 2 separate SNPs, maybe needs to change but how often is this going to be an issue
    return len(h1_only) + len(h2_only)

def sum_base_support(dict1,dict2):
    out_dict = {}
    for key in dict1:
        out_dict[key] = dict1[key] + dict2[key]
    return out_dict

def calc_sd(haplotypes):
    haploCount = 1
    evenDict = {}
    for h in haplotypes:
        #get list of read support vals
        agreeReads = []
        for pos in h:
            reads = h[pos]['base support counts'][h[pos]['assigned allele']]
            agreeReads.append(reads)

        mean = np.mean(agreeReads)
        norm = [x/mean for x in agreeReads]
        evenDict[haploCount] = np.std(norm)
        haploCount += 1
    return evenDict