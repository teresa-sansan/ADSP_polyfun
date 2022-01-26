
## chr9,12,18,22 have missing file in p<5e-8

#cd /gpfs/commons/home/tlin/data/plink_tutorial/height_test
cd /gpfs/commons/home/tlin/data/plink_tutorial/height_test/only_change_PLINK_input

#for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
for i in 0.3 0.2
do
echo write pT_$i.prs
awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' pT_*.$i.profile > pT_$i.prs
done


