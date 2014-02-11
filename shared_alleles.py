#NEANDERTHAL GENOME PROJECT
#Comparative analysis of neanderthal and human subgroup genomes
#Brian Conroy


import random
import numpy


def loadNeandReadD(neandFileName):
    #Returns hash of (chromosome,position):nucleotide values 
    f = open(neandFileName) #file of chromosomes, positions, and alleles
    D ={}   
    while True:
        s = f.readline()
        if s=="":
            break
        chrom,pos,alleles = s.split()
        splitalleles = alleles.split(",")
        D[(chrom,pos)] = tuple(splitalleles)
    return D


def loadModernHumanData(fileName, neandReadD):
    #makes a list of human chromosome/position/allele data for human DNA
    f = open(fileName)
    s = f.readline()
    header = s.split()
    modernHumanData = []
    while True:
        s = f.readline()
        if s == "":
            break
        s = s.rstrip()
        strL= s.split('\t')
        Chr=strL[0]
        pos=strL[1]
        chimp=strL[4]
        genDataL=strL[5:]
        if (Chr,pos) in neandReadD: #so data corresponds to Neanderthal 
            genotypesL = []
            for dipgeno in genDataL:
                a,b=dipgeno.split(",")
                genotypesL.append(a)
                genotypesL.append(b)
            modernHumanData.append((Chr,pos,chimp,tuple(genotypesL)))
    return modernHumanData
                
       
def derAlleleCount(neandReadD,h1SiteDataL,h2SiteDataL):
    "Given a dictionary of neanderthal reads and site data for 2 human"
    "samples, returns a tuple of counts of the number of times h1 and h2 share"
    "derived alleles with neanderthal, respectively."
    count1 = 0
    count2 = 0
    for i in range(0,len(h1SiteDataL)):
        chrom = h1SiteDataL[i][0]
        pos = h1SiteDataL[i][1]
        chimp = h1SiteDataL[i][2]
        randomH1 = random.choice(h1SiteDataL[i][3])
        randomH2 = random.choice(h2SiteDataL[i][3])
        randomN = random.choice(neandReadD[(chrom,pos)])
        if randomH1 != randomH2: #if the human alleles differ
            if randomN != chimp: #if the neanderthal has the derived allele
                if randomH1 == randomN: #if randomH1 shares the derived allele
                    if randomH2 == chimp: #assures there are only 2 allele types
                        count1 += 1
                if randomH2 == randomN: #if randomH2 shares the derived allele
                    if randomH1 == chimp: #same assurance as above
                        count2 += 1                    
    return (count1,count2)
                    
                    
def analysis(): #Takes about 6 minutes to run on my machine. 
    #Compares counts of alleles shared with Neanderthals among Finns,
    #Yorubas, and Chinese 
    D = loadNeandReadD("SLVi33.16.tsv")
    fin = loadModernHumanData("fin-2.10000000-25000000.tsv",D)
    yor = loadModernHumanData("yor-2.10000000-25000000.tsv",D)
    chi = loadModernHumanData("chb-2.10000000-25000000.tsv",D)

    finyorL = []
    for i in range(0,100): #fin - yor

        finCount = derAlleleCount(D,fin,yor)[0]
        yorCount = derAlleleCount(D,fin,yor)[1]
        diff = finCount - yorCount
        finyorL.append(diff)
    meanfinyor = numpy.mean(finyorL)
    stdfinyor = numpy.std(finyorL)

    finchiL = []
    for i in range(0,100): #fin - chi

        finCount = derAlleleCount(D,fin,chi)[0]
        chiCount = derAlleleCount(D,fin,chi)[1]
        diff = finCount - chiCount
        finchiL.append(diff)
    meanfinchi = numpy.mean(finchiL)
    stdfinchi = numpy.std(finchiL)

    yorchiL = []
    for i in range(0,100): #fin - chi

        yorCount = derAlleleCount(D,yor,chi)[0]
        chiCount = derAlleleCount(D,yor,chi)[1]
        diff = yorCount - chiCount
        yorchiL.append(diff)
    meanyorchi = numpy.mean(yorchiL)
    stdyorchi = numpy.std(yorchiL)

    print "For Finnish and Yoruba, mean is"
    print meanfinyor
    print "and stdev is"
    print stdfinyor
    print "For Finnish and Chinese, mean is"
    print meanfinchi
    print "and stdev is"
    print stdfinchi
    print "For Yoruba and Chinese, mean is"
    print meanyorchi
    print"and stdev is"
    print stdyorchi
    

                
        
       
        
        
