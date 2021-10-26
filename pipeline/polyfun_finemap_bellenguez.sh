#!/bin/bash
#SBATCH --job-name=bellenguez_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=28:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_susie/max_snp_1/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/bellenguez_bl.${chr}.snpvar_constrained.gz"
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all.${chr}.snpvar_ridge_constrained.gz"
#sumstat="output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac.${chr}.snpvar_ridge_constrained.gz"


###susie
for i in $FILES
do	
	#echo $i
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	python finemapper.py \
		--sumstats $sumstat \
		--n 487511 \
		--non-func \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 1 \
	  	--allow-missing \
		--out "/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_susie/max_snp_1/finemap_bellenguez_susie.${chr}.$start.$end.gz"

	

done
