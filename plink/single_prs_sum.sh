## This script is for summing the CHR-separated PRS into one PRS value for each individual.
## Because we cauclated thiem by chr to speed up
## The diff between this script and tje prs_pT_sum.py is that
## This scirpt is design for only one obervation (while the other are summing PRS based on different p value thresholds)

#cd /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP
#file_name='updateRSID_interested_SNP_minPIP_chr'
#output_name='updateRSID_interested_SNP_minPIP'

cd /gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta
file_name='max_snp_10_clump_pT_chr'
output_name='new_beta_max_snp_10_clump_pT'

awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $file_name*.profile > $output_name.prs
echo 'finished, write sum up PRS in' ${output_name}.prs
