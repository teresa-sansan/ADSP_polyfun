## This script is for summing the CHR-separated PRS into one PRS value for each individual.
## Because we cauclated thiem by chr to speed up
## The diff between this script and the prs_pT_sum.py is that
## This scirpt is design for only one obervation (while the other are summing PRS based on different p value thresholds)



path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/'
cd $path

# for thres in e-5 0.001 0.005 0.01 0.05 0.1 0.5
# do
#     awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr*.${thres}.profile  > $path/plink_output/prs_${thres}.tsv
# done



#awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr12.chunk{1,2,3}_0.5.prs.tsv > $path/plink_output/test


thres=0.5
echo $path/plink_output/chr12*.${thres}.profile
#awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr12*.${thres}.profile  > $path/plink_output/test