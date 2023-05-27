## This script is for summing the CHR-separated PRS into one PRS value for each individual.
## Because we cauclated thiem by chr to speed up
## The diff between this script and tje prs_pT_sum.py is that
## This scirpt is design for only one obervation (while the other are summing PRS based on different p value thresholds)



path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/'
cd $path

cd /gpfs/commons/home/tlin/output/cT/chr_sep_plink/wightman/new_beta
file_name='max_snp_10_clump_pT_chr'
output_name='new_beta_max_snp_10_clump_pT'

awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' plink_output/$file_name*.profile > $output_name.prs
echo 'finished, write sum up PRS in' ${output_name}.prs

_0.001.prs.tsv

chr9.chunk2_0.001.prs.tsv 
for thres in e-5 0.001 0.005 0.01 0.05 0.1 0.5
do
    awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr*_${thres}.prs.tsv  > $path/plink_output/prs_${thres}.tsv
done