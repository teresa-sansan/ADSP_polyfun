#!/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1> /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/pT/max_10.log 2>&1

#for chr in {2..3}
#do
#echo start chr $chr
#agg_file="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/chr${chr}.aggregrate.all.txt"
#gunzip $agg_file.gz

#awk '{print $2,$10}' $agg_file  > SNP.pvalue

~/plink \
--bfile /gpfs/commons/home/tlin/data/biallelic/${chr}_filt \
--score ${agg_file} 2 4 12 'header' --q-score-range /gpfs/commons/home/tlin/polyfun_script/pipeline/range_list.txt  SNP.pvalue \
--extract /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/clumping/clump_max10_chr${chr}.clumped \
--out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/pT/max_snp_10/pT_max10_merge_chr${chr}

done


##$2 for SNP_ID,$4 for effective allele info, $12 for effect size estimate 



#~/plink --bfile /gpfs/commons/home/tlin/data/biallelic/22_filt --score /gpfs/commons/home/tlin/ou
#t/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/chr22.aggregrate.all.txt 2 4 12 --q-score-range /gpfs/commons/home/tli
#n/polyfun_script/pipeline/range_list.txt SNP.pvalue --extract /gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_co
#nstrained/clumping/clump_max10_chr22.clumped --out pt_max10  
