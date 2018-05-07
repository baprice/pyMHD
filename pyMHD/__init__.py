import subprocess

def check_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which

    return which(name) is not None


def searchKmer(chrom, start, ref, genome, flankSize = 100, trfdir = None):
    """
    Search for microhomolgous kmer around given deletion 
    
    Parameters
    ----------
    chrom : string
        Chromosome location of variant (e.g. 'chr10')
    start : int
        Start position of variant
    ref : string
        Reference allele of variant
    genome : string
        Location of reference genome
    flankSize : int
        Distance around deletion to check for microhomology [100]
    trfdir : string
        Location of TRF (http://tandem.bu.edu/trf/trf.unix.help.html) script. If not set, then microsatellites will not be filtered. 
        
    Returns
    -------
    string
        Microhomolgous kmer.  -1 if deletion is identified as short tandem repeat
        
    """
    
    #Check if samtools is installed
    if not check_tool('samtools'):
        print("ERROR: samtools not found in PATH")
        print("Please install samtools and add to PATH")
        return

    #Get flanking sequence using samtools
    sequence = chrom + ':' + str(int(start)-flankSize) + '-' + str(int(start)+len(ref)+flankSize)
    cmd = ' '.join(['samtools','faidx',genome,sequence])
    up = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    flankseq = ''
    upstream = ''
    for fline in [str(x, 'utf-8') for x in up.stdout.readlines()]:
        if fline[0] == '>':
            continue
        else:
            flankseq = flankseq + fline.strip()
    retvalUp = up.wait()
    upstream = flankseq[:flankSize]
    downstream = flankseq[flankSize+len(ref):]

    #Check for repeats if trfdir specified
    #http://tandem.bu.edu/trf/trf.unix.help.html

    if trfdir is not None:
        trfCmd = 'echo ">seq1\n' + flankseq + '" | ' + trfdir + ' - 2 7 7 80 10 50 7 -h -ngs'
        trf = subprocess.Popen(trfCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if len(trf.stdout.readlines()) > 0:
            return -1
        
    #Construct microhomologous kmer
    #Check upstream for kmer
    kmer = ''
    for i in range(-1,-(len(ref)+1),-1):
        if ref[i] == upstream[i]:
            kmer = kmer + str(ref[i])
        else:
            break
    kmer = kmer[::-1]

    if kmer == '': 
        #Check downstream for kmer
        for i in range(len(ref)):
            if ref[i] == downstream[i]:
                kmer = kmer + str(ref[i])
            else:
                break
        kmer = kmer[::-1]
    
    return kmer



def mafMHD(maf, genome, minDelSize = 5, flankSize = 100, trfdir = None):
    """
    Count MHD in MAF file
    
    Parameters
    ----------
    maf : string
        Location of input MAF file
    genome : string
        Location of reference genome fasta
    minDelSize : int
        Minimum bp deletion size to be considered [5]
    flankSize : int
        Distance around deletion to check for microhomology [100]
    trfdir : string
        Location of TRF (http://tandem.bu.edu/trf/trf.unix.help.html) script. If not set, then microsatellites will not be filtered.
    
    Returns
    -------
    dict
        Dictionary of MHD counts by required size of microhomology
    dict
        File Info: Number of variants, deletions, filtered deletions, short tandem repeats
    """
    

    with open(maf, 'r') as file:
        mhdCount = dict.fromkeys(list(range(0,10)) + ['10+'], 0)
        variantCount = 0
        delCount = 0
        delFiltered = 0
        trfCount = 0
        
        for line in file:
            #Skip line if header
            if line.startswith('#'):
                continue
            
            #Find indicies of required columns
            lineSplit = line.split('\t')
            if 'Hugo_Symbol' in lineSplit:
                variantTypeInd = lineSplit.index('Variant_Type')
                chromInd = lineSplit.index('Chromosome')
                startInd = lineSplit.index('Start_position')
                endInd = lineSplit.index('End_position')
                refInd = lineSplit.index('Reference_Allele')
                continue
            
            variantCount += 1
            try:
                variantType = lineSplit[variantTypeInd]
                chrom = lineSplit[chromInd]
                start = lineSplit[startInd]
                end = lineSplit[endInd]
                ref = lineSplit[refInd]
            except:
                print("ERROR: Required columns not found in MAF file")
                print("Please Check file for [Chromosome, Start_position, End_position, Reference_allele, Variant_Type] columns")
                return
            
            #If variant is not a deletion, skip
            if variantType != 'DEL':
                continue
            else:
                delCount += 1
                             
            #Format chromosome name for samtools
            if 'chr' not in chrom:
                chrom = 'chr' + str(chrom)
                
            #Check deletion size
            if len(ref) < minDelSize:
                delFiltered += 1
                continue
            
            #Check for microhomologous kmer
            kmer = searchKmer(chrom, start, ref, genome, flankSize, trfdir)
            
            #tandem repeat if kmer = -1, skip
            if kmer == -1:
                trfCount += 1
                continue
            
            #Check length of kmer and add to mhdCount
            if len(kmer) < 10:
                for i in range(0, len(kmer)+1):
                    mhdCount[i] += 1
            else:
                mhdCount['10+'] += 1       
            
    sampleInfo = {'Total Variants':variantCount, 'Deletions':delCount, 'minDelFiltered':delFiltered, 'Short Tandem Repeat Filtered':trfCount}
    return mhdCount, sampleInfo


def vcfMHD(vcf, genome, minDelSize = 5, flankSize = 100, trfdir = None):
    """
    Count MHD in VCF file
    
    Parameters
    ----------
    maf : string
        Location of input VCF file
    genome : string
        Location of reference genome fasta
    minDelSize : int
        Minimum bp deletion size to be considered [5]
    flankSize : int
        Distance around deletion to check for microhomology [50]
    trfdir : string
        Location of TRF (http://tandem.bu.edu/trf/trf.unix.help.html) script. If not set, then microsatellites will not be filtered.
    
    Returns
    -------
    dict
        Dictionary of MHD counts by required size of microhomology
    dict
        File Info: Number of variants, deletions, filtered deletions, short tandem repeats
    """
    from cyvcf2 import VCF
    
    #Load in vcf
    vcf = VCF(vcf_file,gts012=True)

    #get FORMAT field list
    format_fields = []
    for h in vcf.header_iter():
        if h.info().get('HeaderType') == 'FORMAT':
            format_fields.append(h.info().get('ID'))

    #get INFO field list
    info_fields = []
    for h in vcf.header_iter():
        if h.info().get('HeaderType') == 'INFO':
            info_fields.append(h.info().get('ID'))
    
    variantCount = 0
    mhdCount = dict.fromkeys(list(range(0,10)) + ['10+'], 0)
    variantCount = 0
    delCount = 0
    delFiltered = 0
    trfCount = 0    
    for variant in vcf:
        variantCount += 1
        altNum = 0
        for alt in variant.ALT:
            altNum += 1
            chrom = str(variant.CHROM)
            start = str(variant.POS)
            ref = str(variant.REF)
                             
            #Format chromosome name for samtools
            if 'chr' not in chrom:
                chrom = 'chr' + str(chrom)
                
            #Check if deletion 
            if len(ref) > len(alt):
                delCount += 1
            else:
                continue
                
            #Check if deletion is large enough
            if len(ref) < minDelSize:
                delFiltered += 1
                continue
                
            #Check for microhomologous kmer
            kmer = searchKmer(chrom, start, ref, genome, flankSize, trfdir)
            
            #tandem repeat if kmer = -1, skip
            if kmer == -1:
                trfCount += 1
                continue

            #Check length of kmer and add to mhdCount
            if len(kmer) < 10:
                for i in range(0, len(kmer)+1):
                    mhdCount[i] += 1
            else:
                mhdCount['10+'] += 1                  
            
    sampleInfo = {'Total Variants':variantCount, 'Deletions':delCount, 'minDelFiltered':delFiltered, 'Short Tandem Repeat Filtered':trfCount}
    return mhdCount, sampleInfo 

