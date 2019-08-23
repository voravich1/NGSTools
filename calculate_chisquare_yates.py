import numpy as np
from scipy.stats import chi2_contingency

treatFile_infect = "/Volumes/4TBSeagateBackupPlus_Drive/opas_bgi_set2/analysis/res_qc80_pass/res_map_virus_bwa_single_summary_strand/5.6_map_sorted_all_htseqCount_strandno.txt"
treatFile_uninfect = "/Volumes/4TBSeagateBackupPlus_Drive/opas_bgi_set2/analysis/res_qc80_pass/res_map_virus_bwa_single_summary_strand/5.5_map_sorted_all_htseqCount_strandno.txt"
untreatFile_infect = "/Volumes/4TBSeagateBackupPlus_Drive/opas_bgi_set2/analysis/res_qc80_pass/res_map_virus_bwa_single_summary_strand/5.8_map_sorted_all_htseqCount_strandno.txt"
untreatFile_uninfect = "/Volumes/4TBSeagateBackupPlus_Drive/opas_bgi_set2/analysis/res_qc80_pass/res_map_virus_bwa_single_summary_strand/5.7_map_sorted_all_htseqCount_strandno.txt"
saveFileName = "/Volumes/4TBSeagateBackupPlus_Drive/opas_bgi_set2/analysis/res_qc80_pass/res_map_virus_bwa_single_summary_strand/rep3_chi_square_yates.txt"

saveFile = open(saveFileName,'w')

treat_infect = np.loadtxt(fname=treatFile_infect,dtype="str")
treat_uninfect = np.loadtxt(fname=treatFile_uninfect,dtype="str")
untreat_infect = np.loadtxt(fname=untreatFile_infect,dtype="str")
untreat_uninfect = np.loadtxt(fname=untreatFile_uninfect,dtype="str")

saveFile.write("POS\tChi_square\tP_Value\tDegree_freedom\tExp_Val_t_f\tExp_Val_t_uf\tExp_Val_ut_f\tExp_Val_ut_uf\n")

for t_f,t_uf,ut_f,ut_uf in zip(treat_infect,treat_uninfect,untreat_infect,untreat_uninfect) :

    t_f_count = int(t_f[1])
    t_uf_count = int(t_uf[1])
    ut_f_count = int(ut_f[1])
    ut_uf_count = int(ut_uf[1])
    try:
        obs = np.array([[t_f_count,t_uf_count],[ut_f_count,ut_uf_count]])
        res = chi2_contingency(obs)
        chi_sq = res[0]
        p_value = res[1]
        degree_freedom = res[2]
        exp_val = res[3]
        #expect_value = "t_f:"+str(exp_val[0,0])+"|t_uf:"+str(exp_val[0,1])+"|ut_f:"+str(exp_val[1,0])+"|ut_uf:"+str(exp_val[1,1])
        #expect_value = str(round(exp_val[0, 0],6)) + "\t" + str(round(exp_val[0, 1],6)) + "\t" + str(round(exp_val[1, 0],6)) + "\t" + str(round(exp_val[1, 1],6))
        expect_value = str(exp_val[0, 0]) + "\t" + str(exp_val[0, 1]) + "\t" + str(exp_val[1, 0]) + "\t" + str(exp_val[1, 1])

        saveFile.write(t_f[0] + "\t" + str(chi_sq) + "\t" + str(p_value) + "\t" + str(degree_freedom) + "\t" + expect_value + "\n")
    except:
        saveFile.write(t_f[0] + "\tNA\tNA\tNA\tNA\n")

saveFile.close()