#! /Users/worawich/miniconda3/envs/NGSTools/bin/python


import sys


snpDBFile = "/Users/worawich/Dropbox/Work/TB_proposal/Document/SNV_lineage_DB_Full.txt"
customAnnFile = "/Users/worawich/Dropbox/Work/TB_proposal/Document/SNV_lineage_DB_Full.vcf"

inputFile = open(snpDBFile,'r')
saveFile = open(customAnnFile,'w')

# Write vcf header
saveFile.write("##fileformat=VCFv4.2\n")
saveFile.write("##reference=/Users/worawich/TBProfiler/ref/MTB-h37rv_asm19595v2-eg18.fa\n")
saveFile.write("##contig=<ID=Chromosome,length=4411532,species=\"Mycobacterium tuberculosis\">\n")
saveFile.write("##INFO=<ID=LIN,Number=1,Type=String,Description=\"Lineage\">\n")
saveFile.write("##INFO=<ID=LOID,Number=1,Type=String,Description=\"Locus ID\">\n")
saveFile.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene name\">\n")
saveFile.write("##INFO=<ID=MUTYPE,Number=1,Type=String,Description=\"Mutation type\">\n")
saveFile.write("##INFO=<ID=GCOOR,Number=1,Type=Integer,Description=\"Gene coordinate\">\n")
saveFile.write("##INFO=<ID=CONUM,Number=1,Type=Integer,Description=\"Codon number\">\n")
saveFile.write("##INFO=<ID=COCHANGE,Number=1,Type=String,Description=\"Codon change\">\n")
saveFile.write("##INFO=<ID=AMCHANGE,Number=1,Type=String,Description=\"Amino acid change\">\n")
saveFile.write("##INFO=<ID=ES,Number=1,Type=String,Description=\"Essentiality\">\n")
saveFile.write("#CHROME\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

headerFlag = True

for line in inputFile:
    if headerFlag == True:
        headerFlag = False
        continue

    info = line.split()
    lineage = info[0]
    pos = info[1]
    gcoor = info[2]
    ref = info[3].split(">")[0]
    alt = info[3].split(">")[1]
    conum = info[4]
    cochange = info[5]
    amchange = info[6]
    loid = info[7]
    gene = info[8]
    mutype = info[9]
    es = info[10]

    infoString = "LIN="+lineage+";"+"GENE="+gene+";"+"LOID="+loid+";"+"MUTYPE="+mutype+";"+"GCOOR="+gcoor+";"+"CONUM="+conum+";"+"COCHANGE="+cochange+";"+"AMCHANGE="+amchange+";"+"ES="+es
    saveFile.write("Chromosome\t"+pos+"\t"+"."+"\t"+ref+"\t"+alt+"\t"+"."+"\t"+"."+"\t"+infoString+"\n")

saveFile.close()
inputFile.close()