import os
import ctk_utils
import numpy as np
import math
from scipy.stats import fisher_exact
from pathlib import Path

input_folder = "C:/Users/vorav/Desktop/opas_bgi_set2/test_post_ctk"
#p = Path("C:/Users/vorav/Desktop/opas_bgi_set2/ctk_post_hg38_n")
p = Path("/Volumes/2TBSeagateBackupPlus/opas_bgi_set2/analysis/new_run/ctk_post_virus_p")

useable_read_pool_inf = p / "pool_inf_virus_p_pool.bed"
useable_read_pool_non = p / "pool_non_virus_p_pool.bed"
useable_read_sm_inf = p / "sm_inf_virus_p_pool.bed"
useable_read_sm_non = p / "sm_non_virus_p_pool.bed"

pool_inf_peak_bed = p / "pool_inf_virus_p_pool_peak_sig.bed"
pool_non_peak_bed = p / "pool_non_virus_p_pool_peak_sig.bed"
#sm_inf_peak_bed = p / "sm_virus_inf_r2_pool_peak_sig.bed"
#sm_non_peak_bed = p / "sm_virus_non_r2_pool_peak_sig.bed"

#union_bedgraph_ctk = "/drives/c/Users/vorav/Desktop/opas_bgi_set2/test_post_ctk/union.bedgraph"
pool_inf_bedgraph = p / "pool_inf_virus_p_pool.bedgraph"
pool_non_bedgraph = p / "pool_non_virus_p_pool.bedgraph"
sm_inf_bedgraph = p / "sm_inf_virus_p_pool.bedgraph"
sm_non_bedgraph = p / "sm_non_virus_p_pool.bedgraph"

# Check peak File is empty or not
if(os.stat(pool_inf_peak_bed).st_size == 0): # Will return true if empty
    print("inf peak file is empty. Can not process")
    exit()
elif(os.stat(pool_non_peak_bed).st_size == 0):
    print("non peak file is empty. Can not process")
    exit()

# Save File name
stat_inf_filename = p / "stat_inf.txt"
stat_non_filename = p / "stat_non.txt"

sig_inf_peak_bedgraph = p / "sig_peak_inf.bedgraph"
sig_non_peak_bedgraph = p / "sig_peak_non.bedgraph"

#useable_read_pool_inf = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_inf_r2_pool.bed"
#useable_read_pool_non = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_non_r2_pool.bed"
#useable_read_sm_inf = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_inf_r2_pool.bed"
#useable_read_sm_non = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_non_r2_pool.bed"

#pool_inf_peak_bed = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_inf_r2_pool_peak_sig.bed"
#pool_non_peak_bed = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_non_r2_pool_peak_sig.bed"
#sm_inf_peak_bed = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_inf_r2_pool_peak_sig.bed"
#sm_non_peak_bed = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_non_r2_pool_peak_sig.bed"

#union_bedgraph_ctk = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/union.bedgraph"
#pool_inf_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_inf_r2_pool.bedgraph"
#pool_non_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/pool_virus_non_r2_pool.bedgraph"
#sm_inf_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_inf_r2_pool.bedgraph"
#sm_non_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sm_virus_non_r2_pool.bedgraph"

#stat_inf_filename = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/stat_inf.txt"
#stat_non_filename = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/stat_non.txt"

#sig_inf_peak_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sig_peak_inf.bedgraph"
#sig_non_peak_bedgraph = "/Users/worawich/Google_Drive/WORK/opas_bgi_set2/ctk_result_for_post_process/sig_peak_non.bedgraph"

stat_inf = open(stat_inf_filename,'w')
stat_non = open(stat_non_filename,'w')
sig_inf_peak = open(sig_inf_peak_bedgraph,'w')
sig_non_peak = open(sig_non_peak_bedgraph,'w')

read_pool_inf = np.loadtxt(fname=useable_read_pool_inf,dtype="str")
total_read_pool_inf = read_pool_inf.shape[0]

read_pool_non = np.loadtxt(fname=useable_read_pool_non,dtype="str")
total_read_pool_non = read_pool_non.shape[0]

read_sm_inf = np.loadtxt(fname=useable_read_sm_inf,dtype="str")
total_read_sm_inf = read_sm_inf.shape[0]

read_sm_non = np.loadtxt(fname=useable_read_sm_non,dtype="str")
total_read_sm_non = read_sm_non.shape[0]

#pool_inf_peak = np.loadtxt(fname=pool_inf_peak_bed,dtype="str")
#pool_non_peak = np.loadtxt(fname=pool_non_peak_bed,dtype="str")

#union_bedgraph = np.loadtxt(fname=union_bedgraph_ctk,dtype="str")

union_count_dict = {}
pool_inf_count_dict = ctk_utils.createCountDictFromBedgraph(pool_inf_bedgraph)
pool_non_count_dict = ctk_utils.createCountDictFromBedgraph(pool_non_bedgraph)
sm_inf_count_dict = ctk_utils.createCountDictFromBedgraph(sm_inf_bedgraph)
sm_non_count_dict = ctk_utils.createCountDictFromBedgraph(sm_non_bedgraph)


#Convert to RPM bedgraph
ctk_utils.convert2rpmBedgraph(pool_inf_bedgraph,total_read_pool_inf)
ctk_utils.convert2rpmBedgraph(pool_non_bedgraph,total_read_pool_non)
ctk_utils.convert2rpmBedgraph(sm_inf_bedgraph,total_read_sm_inf)
ctk_utils.convert2rpmBedgraph(sm_non_bedgraph,total_read_sm_non)

# Loop create universal dict count. Value column represent inf non_inf sm_inf sm_non_inf
#first_flag = True
#for union_data in union_bedgraph:

    #if(first_flag == True):

     #   union_count_dict['header'] = [union_data[3],union_data[4],union_data[5],union_data[6]]
      #  first_flag = False
    #else:
        #start = int(union_data[1])
        #stop = int(union_data[2])

        #for i in range(start,stop+1):
            #union_count_dict[i]=[union_data[3],union_data[4],union_data[5],union_data[6]]

# Prerequisit data
pool_inf_pm_factor = total_read_pool_inf/1000000
pool_non_pm_factor = total_read_pool_non/1000000
sm_inf_pm_factor = total_read_sm_inf/1000000
sm_non_pm_factor = total_read_sm_non/1000000
sm_inf_count = 0
sm_non_count = 0

# Loop inf peak. Fisher exact test + log2 fold change
stat_inf.write("ID\tPOS_START\tPOS_STOP\tlog2FoldChange\tOdds_ratio\tP_Value\n")
sig_inf_peak.write("track type=bedGraph name=\"Significant peak (infect)\"\n")

with open(pool_inf_peak_bed,"r") as pool_inf_peak:

    for line in pool_inf_peak:
        if(line == "\n"):
            break

        peak_data = line.split("\t")
        chromosome = peak_data[0]
        start = int(peak_data[1])
        stop = int(peak_data[2])
        peak_id = peak_data[3]
        treat_inf_peak = int(peak_data[4])

        treat_inf_other = total_read_pool_inf - treat_inf_peak
        sm_inf_count = 0
        for pos in range(start,stop):

            try:
                count = sm_inf_count_dict[pos]
                sm_inf_count = sm_inf_count + count
            except:
                count=0

        sm_inf_peak = round(sm_inf_count/(stop-start))
        sm_inf_other = total_read_sm_inf - sm_inf_peak

        # Fisher exact test
        input_array_table = np.array([[treat_inf_peak, treat_inf_other], [sm_inf_peak, sm_inf_other]])
        res = fisher_exact(input_array_table)
        oddsratio = res[0]
        p_value = res[1]
        #degree_freedom = res[2]
        #exp_val = res[3]
        # expect_value = "t_f:"+str(exp_val[0,0])+"|t_uf:"+str(exp_val[0,1])+"|ut_f:"+str(exp_val[1,0])+"|ut_uf:"+str(exp_val[1,1])
        # expect_value = str(round(exp_val[0, 0],6)) + "\t" + str(round(exp_val[0, 1],6)) + "\t" + str(round(exp_val[1, 0],6)) + "\t" + str(round(exp_val[1, 1],6))
        #expect_value = str(exp_val[0, 0]) + "\t" + str(exp_val[0, 1]) + "\t" + str(exp_val[1, 0]) + "\t" + str(exp_val[1, 1])

        # log2 fold change
        rpm_treat_inf_peak = treat_inf_peak/pool_inf_pm_factor
        rpm_sm_inf_peak = sm_inf_peak/sm_inf_pm_factor
        # check to prevent zero divide
        if(rpm_sm_inf_peak == 0):
            rpm_sm_inf_peak = 1
        ratio = rpm_treat_inf_peak/rpm_sm_inf_peak
        l2fold = round(math.log(ratio,2))

        stat_inf.write(peak_id + "\t" + str(start) + "\t" + str(stop) + "\t"+ str(l2fold) + "\t" + str(oddsratio) + "\t" + str(p_value) + "\n")

        if(p_value<0.05):
            sig_inf_peak.write(chromosome + "\t" + str(start) + "\t" + str(stop) + "\t" + str(round((rpm_treat_inf_peak)/100000,2)) + "\n")  # we scale peak count down to 10^5 scale

stat_inf.close()
sig_inf_peak.close()

########################################################

# Loop non peak. Fisher exact test + log2 fold change
stat_non.write("ID\tPOS_START\tPOS_STOP\tlog2FoldChange\tOdds_ratio\tP_Value\n")
sig_non_peak.write("track type=bedGraph name=\"Significant peak (non-infect)\"\n")

with open(pool_non_peak_bed,"r") as pool_non_peak:
    for line in pool_non_peak:
        peak_data = line.split("\t")
        if (line == "\n"):
            break

        peak_data = line.split("\t")
        chromosome = peak_data[0]
        start = int(peak_data[1])
        stop = int(peak_data[2])
        peak_id = peak_data[3]
        treat_non_peak = int(peak_data[4])

        treat_non_other = total_read_pool_non - treat_non_peak
        sm_non_count = 0
        for pos in range(start,stop):

            try:
                count = sm_non_count_dict[pos]
                sm_non_count = sm_non_count + count
            except:
                count=0

        sm_non_peak = round(sm_non_count/(stop-start))
        sm_non_other = total_read_sm_non - sm_non_peak

        # Fisher exact test
        input_array_table = np.array([[treat_non_peak, treat_non_other], [sm_non_peak, sm_non_other]])
        res = fisher_exact(input_array_table)
        oddsratio = res[0]
        p_value = res[1]
        #degree_freedom = res[2]
        #exp_val = res[3]
        # expect_value = "t_f:"+str(exp_val[0,0])+"|t_uf:"+str(exp_val[0,1])+"|ut_f:"+str(exp_val[1,0])+"|ut_uf:"+str(exp_val[1,1])
        # expect_value = str(round(exp_val[0, 0],6)) + "\t" + str(round(exp_val[0, 1],6)) + "\t" + str(round(exp_val[1, 0],6)) + "\t" + str(round(exp_val[1, 1],6))
        #expect_value = str(exp_val[0, 0]) + "\t" + str(exp_val[0, 1]) + "\t" + str(exp_val[1, 0]) + "\t" + str(exp_val[1, 1])

        # log2 fold change
        rpm_treat_non_peak = treat_non_peak/pool_non_pm_factor
        rpm_sm_non_peak = sm_non_peak/sm_non_pm_factor
        # check to prevent zero divide
        if(rpm_sm_non_peak == 0):
            rpm_sm_non_peak = 1
        ratio = rpm_treat_non_peak/rpm_sm_non_peak
        l2fold = round(math.log(ratio,2))

        stat_non.write(peak_id + "\t" + str(start) + "\t" + str(stop) + "\t"+ str(l2fold) + "\t" + str(oddsratio) + "\t" + str(p_value) + "\n")

        if(p_value<0.05):
            sig_non_peak.write(chromosome + "\t" + str(start) + "\t" + str(stop) + "\t" + str(round((rpm_treat_non_peak)/100000,2)) + "\n")  # we scale peak count down to 10^5 scale

stat_non.close()
sig_non_peak.close()

########################################################


