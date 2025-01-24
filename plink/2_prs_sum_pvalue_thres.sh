dir_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/clumping_pt_oct/prs_middlefile/'
file_name='prs_plink_pT'
cd $dir_path

## summing up per-chr pRS to genomewide PRS
for i in e-6 e-5 e-4 0.001 0.01 0.1 
do
	echo write pT_$i.prs
	awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' ${file_name}_*.$i.profile > pT_${i}.prs
done

# source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
# conda activate polyfun

## Merge with phenotype file
echo 'merge with phenotype file...'
python /gpfs/commons/home/tlin/script/plink/merge_pvalue_diagnosis.py $dir_path /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/c+pt/c_pt.tsv

echo 'done'