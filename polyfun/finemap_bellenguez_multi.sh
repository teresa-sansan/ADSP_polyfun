#!/bin/bash
#SBATCH --job-name=bellenguez_finemap_multi
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=30:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/bellenguez_bl.${chr}.snpvar_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all.${chr}.snpvar_ridge_constrained.gz"
#sumstat="output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac.${chr}.snpvar_ridge_constrained.gz"
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.${chr}.snpvar_constrained.gz"

for i in $FILES
do	
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	python finemapper.py \
		--ld $i \
		--sumstats $sumstat \
		--n 487511 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal $max_num_snp \
	  	--allow-missing \
		--out "/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_${max_num_snp}/finemap_bellenguez.${chr}.$start.$end.gz"

	

done
