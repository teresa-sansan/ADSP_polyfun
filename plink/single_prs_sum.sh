#cd /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP
#file_name='updateRSID_interested_SNP_minPIP_chr'
#output_name='updateRSID_interested_SNP_minPIP'

cd /gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta
file_name='max_snp_10_clump_pT_chr'
output_name='new_beta_max_snp_10_clump_pT'

awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $file_name*.profile > $output_name.prs
echo 'finished, write sum up PRS in' ${output_name}.prs
