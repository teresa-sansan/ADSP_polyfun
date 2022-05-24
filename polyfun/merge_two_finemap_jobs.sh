ori='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10'
new='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/aggregate_all_position'

converge='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/overlap_snp/converge_position.tsv'

#sort $new.tsv > $new.sort.tsv
touch $converge
for i in {1..22}
do
echo start chr $i ...
#cat $ori/chr$i.aggregrate.all.txt| cut -f 1-5| sort > $ori/sort_chr${i}.aggregrate.all.txt
comm -1 $ori/sort_chr${i}.aggregrate.all.txt $new.sort.tsv -3 >> /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge/overlap_snp/uniq_chr${i}.tsv
#comm -1 $ori/sort_chr${i}.aggregrate.all.txt $new.sort.tsv -3 >> $converge

done
