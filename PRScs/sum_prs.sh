'''
This script is for summing the CHR-separated PRS into one PRS value for each individual.
Because we cauclated thiem by chunking chr to speed up
'''


#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/'
#path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez'
path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_ibd_adsp_fixed/bellenguez'
cd $path

for thres in  0.01 0.05 0.1 0.5
do
    echo running $thres
    for chr in {1..22}; do tr -s ' ' < chr${chr}.qc.${thres}.profile | cut -d ' ' -f 2-7 > prscs_chr${chr}_prs.${thres}.tsv; done
    awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/prscs_chr*_prs.${thres}.tsv > $path/prs_${thres}.tsv
    #awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr*.${thres}.profile  > $path/plink_output/prs_${thres}.tsv
done

echo "write everything in $path"
#awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $path/plink_output/chr12*.${thres}.profile  > $path/plink_output/test

