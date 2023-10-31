def parseVCFBases(vcfFile):
    codeToBase = {}
    with open(vcfFile) as file:
        for line in file:
            if line.startswith('#'):
                continue
            else:
                line = line.strip('\n').split('\t')
                pos, ref, alt = line[1], line[3], line[4]
                pos = int(pos)
                codeToBase[pos] = {}
                codeToBase[pos]['0'] = ref
                if ',' in alt:
                    #multiallelic SNP
                    codeToBase[pos]['1'], codeToBase[pos]['2'] = alt.split(',')
                else:
                    codeToBase[pos]['1'] = alt 
    return codeToBase

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

                        
def parseFlopp(flopOutputFile,vcfFile):
    #ughhhh, going to have to parse the vcf as well to get bases
    vDict = parseVCFBases(vcfFile)
    with open(flopOutputFile) as file:
        file.readline() #header, gives ref mapped to
        haps = []
        for line in file:
            if line.startswith('1:'):
                #first line
                line = [x for x in line.strip('\n').split('\t') if x != '']
                k = int((len(line) - 1)/2) #number of haplotypes
                pos = int(line[0].split(':')[1])
                hapCols = line[1:k+1]
                for h in hapCols:
                    #set up haplotype storage
                    haps.append([(pos,vDict[pos][h])])
            elif line == '*****\n':
                #EOF
                continue
            else:
                line = line.strip('\n').split('\t')
                pos = int(line[0].split(':')[1])
                hapCols = line[1:k+1]
                for i in range(k):
                    haps[i].append((pos,vDict[pos][hapCols[i]]))
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

if __name__ == '__main__':
    from sys import argv
    # floppOut = parseFlopp('flopp/floppResultsTriple.txt','flopp/freeBayesTripleForFlopp.vcf')
    # floppOut = parseFlopp(argv[1],argv[2])
    vcfOut = parseVCFCalls(argv[1:-3])
    print(vcfOut)
    makeHaploFasta(vcfOut,argv[-3],argv[-2],argv[-1])
