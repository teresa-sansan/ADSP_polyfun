dir_path='/gpfs/commons/home/tlin/output/wightman/prscs/all_anno'
cd $dir_path
for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
do
cat prscs.${i}.profile |tr -s ' '| cut -d ' ' -f 2-7 > prscs_${i}_prs.tsv
done