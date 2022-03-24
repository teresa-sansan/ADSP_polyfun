#!/bin/bash
#SBATCH --job-name=finemap_susie
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=25:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_susie/%x%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr${chr}*"
sumstat='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.${chr}.snpvar_constrained.gz'

for i in $FILES
do	
	filename=$(echo $i | cut -d'/' -f 10 |cut -d'.' -f1)
	start=$(echo $filename| cut -d'_' -f 2)
	end=$(echo $filename| cut -d'_' -f 3)
	
	if [ $max_num_snp -eq 1 ] ; 
	then
		echo /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_susie/max_snp_${max_num_snp}/$chr
	else	
		echo /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_susie/max_snp_${max_num_snp}/$chr
	fi
done
