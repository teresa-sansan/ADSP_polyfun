
## chr9,12,18,22 have missing file in p<5e-8

cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT

for i in e-5 0.001 0.05 0.01 0.1 0.2 0.5
do
echo write pT_check_$i.prs
awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' pT_check_chr*.$i.profile > pT_check_$i.prs
done


