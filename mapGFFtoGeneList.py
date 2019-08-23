#! /Users/worawich/miniconda3/envs/NGSTools/bin/python

import sys

gene_file = sys.argv[1]
gff_file = sys.argv[2]
map_result_file = sys.argv[3]

#gene_file = "/Volumes/4TB_WD/Opas_NGS_results/target_gene.txt"
#gff_file = "/Volumes/4TB_WD/Opas_NGS_results/clean_res_starIndex100_full/binding_sites_report_hg38/CHKSEI85218110013Aligned_sorted_markDup_binding_sites_report_hg38.txt"
#map_result_file = "/Volumes/4TB_WD/Opas_NGS_results/clean_res_starIndex100_full/binding_sites_report_hg38/CHKSEI85218110013_mapTarget.txt"

targetGeneDict = {}

# open gene list file (2 column tab delimit file first is gene symbol, second is ensemble ID)
file = open(gene_file,'r')
for line in file:
    line_noend = line.strip()
    info = line_noend.split("\t")
    gene_symbol = info[0]
    ensemble_ID = info[1]

    targetGeneDict.update({ensemble_ID: gene_symbol})

file.close()

# open gff file (target to gene ensembleID)
savefile = open(map_result_file,'w')
file = open(gff_file,'r')
for line in file:
    line_noend = line.strip()

    gff_info = line_noend.split("\t")
    attribute = gff_info[8]
    attribute_info = attribute.split(";")
    id = attribute_info[0]
    id_info = id.split("=")[1]
    ensemble_ID = id_info.split(":")[1]

    if ensemble_ID in targetGeneDict:
        savefile.write(line_noend)
        savefile.write("\n")

file.close()
savefile.close()