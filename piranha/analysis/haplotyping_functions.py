import pysam
import subprocess
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
    
