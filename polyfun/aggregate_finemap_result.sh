#!/bin/bash
#SBATCH --job-name=aggregate_finemap
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/%x_%j.log

## Note:
## There are two parts of this script. 
## (1) Aggregate the result from polyfun to aggregate the overlapping regions.  (3 MB window size, 1 MB overlap, iteration = 100)
## (2) Aggregate the result from not converging regions (Those failed to converge with IBSS, 1 MB window size, 0.5 MB overlap. iteration = 1000)


cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

##bellenguez_fixed_0224
if true;then
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sep22_new_bellenguez_et_al_2021_hg37_ldsc.munged.parquet'
sumstat_name='bellenguez'
#path='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations'
#path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/'
path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
fi

##wightman
if false;then
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.munged.parquet'
sumstat_name='wightman'
#prefix='all_anno'
path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'
#path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/'
#path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/all_anno/finemap'
#path='/gpfs/commons/home/tlin/output/wightman/finemap'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224/susie/finemap_fixed_assertion_susie_iter'
fi

##kunkle
if false;then
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Kunkle_et_al_2019_hg37_ldsc.tsv.gz'
sumstat_name='kunkle'
path='/gpfs/commons/home/tlin/output/kunkle/new_anno/'
fi

## jansen
if false;then
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Jansen_et_al_2019_hg37_ldsc.tsv.gz'
path='/gpfs/commons/home/tlin/output/jansen'
#path='/gpfs/commons/home/tlin/output/jansen/susie/'
#prefix="jansen"
#path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations'
#path='/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations'
fi

#anno='susie'
## set true if you want to calculate normal PolyFun finemap result
##10-22 in max_5
if true; then
for chr in {1..22}
do
	echo 'start aggregating chr' $chr
	python aggregate_finemapper_results_modified.py \
		--out-prefix $path/$anno/finemap/max_snp_${max_snp}/$sumstat_name \
		--sumstats $sumstat\
		--out $path/$anno/finemap/max_snp_${max_snp}/chr${chr}.aggregate.all.txt.gz \
		--allow-missing \
       	--chr $chr

done
fi

## Set true if you want to aggregate the non-converging regions.  
## this part is trying to aggregate the regions that fails to converge at the first time.  

#prefix_converge='finemap_max_snp_3'                                                
## bellenguez_fixed_convergence
#path_converge='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/bl_brain_atac/max_snp_10'
#path_converge='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/all_anno/finemap/max_snp_10/'
## wightman fixed convergence ##070622
#path_converge='/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/try_rescue_not_converge'

#anno='all_anno'
anno='susie'
if false; then
python aggregate_finemapper_results_fixed_convergence.py \
                --out-prefix $path/$anno/max_snp_10/try_rescue_not_converge/$prefix_converge \
                --sumstats $jansen \
                --regions-file $path/$anno/max_snp_10/run_IBSS_not_converge_list.txt \
                --allow-missing \
                --out $path/$anno/max_snp_10/try_rescue_not_converge/aggregate_rescue.all.txt.gz
fi
