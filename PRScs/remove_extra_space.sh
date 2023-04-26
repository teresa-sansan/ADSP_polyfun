dir_path='/gpfs/commons/home/tlin/output/wightman/prscs'

for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
do
cat prs.${i}.profile |tr -s ' '| cut -d ' ' -f 2-7 > p_${i}_prs.tsv
done