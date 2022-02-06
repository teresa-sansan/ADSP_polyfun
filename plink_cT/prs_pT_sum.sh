
## chr9,12,18,22 have missing file in p<5e-8

#cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/pT
#cd /gpfs/commons/home/tlin/output/cT/kunkle

cd /gpfs/commons/home/tlin/output/cT/bellenguez

for i in e-5 0.001 0.005 0.01 0.05 0.1 0.5
do
echo write pT_$i.prs
awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' bellenguez_pT_chr*.$i.profile > pT_$i.prs
done


