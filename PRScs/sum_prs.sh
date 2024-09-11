
note=$(cat << EOF
This script is for summing the CHR-separated PRS into one PRS value for each individual.
Because we cauclated thiem by chunking chr to speed up
EOF
)

#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/prscs/'
#path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez'
#path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_ibd_adsp_fixed/bellenguez'
#path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez_rerun_0909'
path='/gpfs/commons/home/tlin/output/prs/PRSCS/36k_adsp_ld_panel/bellenguez/'
savename='bellenguez_adsp_ld'

cd $path/prs_middlefile/
for thres in e-6 e-5 e-4 0.001 0.01 0.1
do

    echo running $thres
    for chr in {1..22}; do tr -s ' ' < chr${chr}.qc.${thres}.profile | cut -d ' ' -f 2-7 > prscs_chr${chr}_prs.${thres}.tsv; done
    awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' prscs_chr*_prs.${thres}.tsv > prs_${thres}.tsv
done

cp prs_${thres}.tsv $path/prs_${thres}.tsv

echo "merging prs with phenotype file"
python /gpfs/commons/home/tlin/script/PRScs/merge_prs.py $path $savename

# 0.01 0.1 