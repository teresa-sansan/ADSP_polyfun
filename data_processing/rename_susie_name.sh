cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_susie/max_snp_10
for i in *
do 
name=$i
chr=$(echo $name| cut -d '_' -f 3|sed s/"chr"//)
rename=$(echo $name| cut -d '_' -f 4)
filename=bellenguez_susie.$chr.$rename
echo $filename
mv $i $filename 
done   
