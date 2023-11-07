import pysam
import os
import subprocess
import csv
from scipy.stats import chisquare

from piranha.utils.log_colours import green,cyan,yellow
from piranha.utils.config import *

#get_bam_partition and write_contig_headers_vcf based on scripts from flopp 
def get_bam_partition(read_part_file,bam_file,pref_nam):
    # read_part_file = sys.argv[1]
    # bam_file = sys.argv[2]
    # pref_nam = sys.argv[3]
    bam = pysam.AlignmentFile(bam_file)

    read_part = []

    count_i = -1 
    for line in open(read_part_file,'r'):
        if '#' in line:
            read_part.append(set())
            count_i += 1
        else:
            read_part[count_i].add(line.split()[0])

    ploidy = len(read_part)

    obam_files = []
    file_names = []

    for i in range(ploidy):
        file_name = pref_nam+str(i)+".bam"
        file_names.append(file_name)
        obam_files.append(pysam.AlignmentFile(file_name, "wb", template=bam))

    file_names.append(pref_nam+"-not_mapped.bam")
    obam_files.append(pysam.AlignmentFile(pref_nam+"-not_mapped.bam", "wb", template=bam))

    for b in bam.fetch(until_eof=True):
        not_frag = True;
        for i in range(ploidy):
            qnames = read_part[i]
            if b.query_name in qnames:
                obam_files[i].write(b)
                not_frag=False
        if not_frag:
            obam_files[-1].write(b);


    for obam in obam_files:
        obam.close()

    for file_name in file_names:
        subprocess.run("samtools index " + file_name, shell=True, check=True)

    bam.close()

def write_contig_headers_vcf(vcf_file,outVCF):
    refs = set()
    for line in open(vcf_file,'r'):
        if line == "":
            continue
        if line[0] == '#':
            continue
        ref_chrom = line.split()[0]
        # print(ref_chrom)
        refs.add(ref_chrom)

    refs = list(refs)
    refs.sort()

    new_vcf =  open(outVCF,'w')
    count = 0
    for line in open(vcf_file,'r'):
        if count != 2:
            new_vcf.write(line)
        else:
            for ref in refs:
                new_vcf.write("##contig=<ID="+ref+">\n")
            new_vcf.write(line)
        count += 1

def parseVCFCalls(vcfFileList):
    #dict form {hapCount:[(pos,base)]}
    haps = []
    i=-1
    for v in vcfFileList:
        i += 1
        haps.append([])
        with open(v) as file:
            for line in file:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip('\n').split('\t')
                    pos, ref, alt = line[1], line[3], line[4]
                    
                    if ',' in alt:
                        pos = int(pos)
                        alts = alt.split(',')
                        for a in alts:
                            haps[i].append((pos,a))

                    elif len(ref) != len(alt):
                        #deletion, skip
                        continue

                    elif len(ref) > 1:
                        #MNP with no del, split
                        haps[i].append((pos,alt[0]))
                        haps[i].append((str(int(pos)+1),alt[1]))

                    else:
                        #just a normal SNP
                        pos = int(pos)
                        haps[i].append((pos,alt))
    return haps

def makeHaploFasta(phasingInfo,fastaOutName,refSeqName,refSeqFile):
    seq = ''
    with open(refSeqFile) as refFile, open(fastaOutName,'w') as out:
        # refHead = refFile.readline()
        refFile.readline()
        refHead = f'>{refSeqName}\n'
        out.write(refHead)
        for line in refFile:
            out.write(line)
            line = line.strip('\n')
            seq += line
        for i in range(len(phasingInfo)):
            haploName = 'Haplotype_' + str(i+1)
            muts = phasingInfo[i]
            inSeq = seq
            for m in muts:
                indToMut = int(m[0])-1
                mutTo = m[1]
                inSeq = inSeq[:indToMut] + mutTo + inSeq[indToMut+1:]
            out.write('>'+haploName+'\n')
            out.write(inSeq+'\n')
            inSeq = seq
    
def parseVCF(vcfFile):
    vcfInfo = {}
    with open(vcfFile) as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                line = line.strip('\n').split('\t')
                pos, ref, alt, info = line[1], line[3], line[4], line[7]
                pos = int(pos)
                infoDict = dict([x.split('=') for x in info.split(';')])
                vcfInfo[pos] = {}
                vcfInfo[pos]['0'] = ref
                if ',' in alt:
                    #multiallelic SNP
                    vcfInfo[pos]['1'], vcfInfo[pos]['2'] = alt.split(',')
                else:
                    vcfInfo[pos]['1'] = alt 
                vcfInfo[pos]['readCounts'] = {}
                refReads = infoDict['RO']
                altReads = infoDict['AO']
                #float because freebayes can have 'fractional' observations
                vcfInfo[pos]['readCounts']['0'] = float(refReads)
                if ',' in altReads:
                    counts = [float(x) for x in altReads.split(',')]
                    vcfInfo[pos]['readCounts']['1'] = counts[0]
                    vcfInfo[pos]['readCounts']['2'] = counts[1]
                else:
                    vcfInfo[pos]['readCounts']['1'] = float(altReads)
                vcfInfo[pos]['frequencies'] = {}
                for allele in vcfInfo[pos]['readCounts']:
                    vcfInfo[pos]['frequencies'][allele] = vcfInfo[pos]['readCounts'][allele]/sum(vcfInfo[pos]['readCounts'].values())
                #frequencies list allele freqs in 0,1,2 order

    return vcfInfo

def parseFloppReads(floppPath):
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
                hapCols = line[1:k+1]
                supportCols = line[k+1:]
                
                for i in range(k):
                    haps.append({})
                    outDict = {}
                    outDict['assigned allele'] = hapCols[i]
                    #something has to change here, want reads for all possible alleles
                    outDict['total assigned'], outDict['base support counts'] = parseSupportCol(supportCols[i])
                    haps[i][pos] = outDict

            elif line == '*****\n':
                #EOF
                continue

            else:
                line = line.strip('\n').split('\t')
                pos = line[0].split(':')[1]
                pos=int(pos)
                hapCols = line[1:k+1]
                supportCols = line[k+1:]

                for i in range(k):
                    outDict = {}
                    outDict['assigned allele'] = hapCols[i]
                    outDict['total assigned'], outDict['base support counts'] = parseSupportCol(supportCols[i])
                    haps[i][pos] = outDict
    return haps


def parseSupportCol(col):
    baseSplit = col.split('|')
    totalReads = 0
    baseSupport = {}
    for base in baseSplit:
        baseVal,supportReads = base.split(':')
        # if baseVal == selBaseVal:
        #     baseSupport = int(supportReads)
        baseSupport[baseVal] = int(supportReads)
        totalReads += int(supportReads)
    return (totalReads,baseSupport)

def checkSiteSig(baseSupport,alleleFreqs):
    #check if allele freqs are consistent with random sample from overall reads for this site
    #need to know these are the same order
    #reads might have 100% one allele, need to insert zero in this case, then sort to make sure order correct
    for key in alleleFreqs:
        if key not in baseSupport:
            baseSupport[key] = 0
    readCounts = [baseSupport[k] for k in sorted(baseSupport)]
    exp = [alleleFreqs[alleleFreq]*sum(readCounts) for alleleFreq in alleleFreqs]
    try:
        pval = chisquare(f_obs=readCounts,f_exp=exp).pvalue
    except Exception as e:
        print('Observed values for read counts',readCounts)
        print('Expected values for read counts',exp)
        print('Allele frequencies',alleleFreqs)
        raise ValueError(f'chi squared calculation failed with error, {e}, see above values')
    if pval < 0.01:
        return True
    else:
        return False

def getHaploDist(h1,h2,parsedVCF):
    #returns the number of reliable SNPs separating 2 haplotypes
    sigPos1, sigPos2 = [],[]
    # for position in each haplo, check significance (i.e. are the base proportions differnt from a random sample of total read set), add to list if is
    #need a lookup from vcf parsing to get allelefreqs
    for pos in h1:
        if checkSiteSig(h1[pos]['base support counts'],parsedVCF[pos]['frequencies']) and h1[pos]['assigned allele'] != '0':
            #need to append both position and allele to list
            sigPos1.append((pos,h1[pos]['assigned allele']))
    for pos in h2:
        if checkSiteSig(h2[pos]['base support counts'],parsedVCF[pos]['frequencies']) and h2[pos]['assigned allele'] != '0':
            sigPos2.append((pos,h2[pos]['assigned allele']))
    #how many of the reliable SNPs are only in one
    #set arith?
    h1Only = set(sigPos1) - set(sigPos2)
    h2Only = set(sigPos2) - set(sigPos1)
    #this counts the same site with different alt alleles as 2 separate SNPs, maybe needs to change but how often is this going to be an issue
    return len(h1Only) + len(h2Only)

def sumBaseSupport(dict1,dict2):
    outDict = {}
    for key in dict1:
        outDict[key] = dict1[key] + dict2[key]
    return outDict


def collapseClose(parsedFlopp,distanceCutoff,vcfFile):
    haploPos = 0
    collapsed = []
    output = []
    vcf = parseVCF(vcfFile)
    for h in parsedFlopp:
        collapseIntoCurrent = [haploPos]
        if haploPos in collapsed:
            #don't consider haplotypes already merged into another
            haploPos += 1
            continue
        newH = {}
        foundOthers = False
        
        #for each haplo, check if any others are identical/within distance cutoff, collapse in and combine reads if so, remove collapsed haplo from further consideration
        #TODO add in merge for haplos not sig different from each other
        #check multiple merges will work correctly
        #I'm so dumb, literally all positions are parsed now so it's not an issue for the collapsing
        haploPos2 = haploPos +1
        for h2 in parsedFlopp[haploPos2:]:
            distance = getHaploDist(h,h2,vcf)
            if distance <= distanceCutoff and newH == {}:
                #collapse values and create new haplo 
                #TODO now needs logic for which allele we accept as correct - majority rules is consistent with the rest of the flopp output
                for pos in h:
                    #combine
                    newH[pos] = {}
                    newH[pos]['total assigned'] = h[pos]['total assigned'] + h2[pos]['total assigned']
                    newH[pos]['base support counts'] = sumBaseSupport(h[pos]['base support counts'],h2[pos]['base support counts'])
                    #check if allele assignment is different now that reads are collapsed
                    newH[pos]['assigned allele'] = max(newH[pos]['base support counts'],key=lambda x:newH[pos]['base support counts'][x])
                    
            
                collapsed.append(haploPos2)
                collapseIntoCurrent.append(haploPos2)
                foundOthers = True

            elif distance <= distanceCutoff:
                for pos in h:
                    newH[pos]['total assigned'] = newH[pos]['total assigned'] + h2[pos]['total assigned']
                    newH[pos]['base support counts'] = sumBaseSupport(newH[pos]['base support counts'],h2[pos]['base support counts'])
                    #check if allele assignment is different now that reads are collapsed
                    newH[pos]['assigned allele'] = max(newH[pos]['base support counts'],key=lambda x:newH[pos]['base support counts'][x])
                        
                collapsed.append(haploPos2)
            haploPos2 += 1

        if not foundOthers:
            newH = h
        haploPos += 1
        output.append(collapseIntoCurrent)
        # outputHaplos.append(newH)
    #change to return which haplos are to be merged 
    return output