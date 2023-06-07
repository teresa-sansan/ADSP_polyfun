
#path='/gpfs/commons/home/tlin/output/wightman/prscs/all_anno'
#path='/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/PIP_not0/'
path='/gpfs/commons/home/tlin/output/wightman/prscs/all_except_enformer/'
cd $path
for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
do
cat prs.${i}.profile |tr -s ' '| cut -d ' ' -f 2-7 > prscs_${i}_prs.tsv
done