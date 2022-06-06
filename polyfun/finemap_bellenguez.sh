#!/bin/bash
#SBATCH --job-name=bellenguez_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=200G
#SBATCH --time=58:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo chr $chr
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/bellenguez_bl.${chr}.snpvar_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all.${chr}.snpvar_ridge_constrained.gz"
#sumstat="output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac.${chr}.snpvar_ridge_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.${chr}.snpvar_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_qc/bellenguez.${chr}.snpvar_constrained.gz"
#sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.${chr}.snpvar_constrained.gz"

sumstat_path="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/"
for i in $FILES
do	
	#echo $i
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	if [ $max_num_snp -eq 1 ]
	then
	python finemapper.py \
		--sumstats $sumstat_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
		--n 487511 \
	  	--chr ${chr} --start $start --end $end \
	  	--method susie \
     	  	--max-num-causal 1 \
	  	--allow-missing \
		--out "/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/$anno/max_snp_1/finemap_bellenguez.${chr}.$start.$end.gz"
	else
	python finemapper.py \
		--ld $i \
                --sumstats $sumstat_path/$anno/${anno}.${chr}.snpvar_constrained.gz \
                --n 487511 \
                --chr ${chr} --start $start --end $end \
                --method susie \
                --max-num-causal $max_num_snp \
                --allow-missing \
                --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_annotations/$anno/max_snp_${max_num_snp}/finemap_bellenguez.${chr}.$start.$end.gz
	fi
done

