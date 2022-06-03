#!/bin/bash
#SBATCH --job-name=aggregate_result
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap/%x_%j.log

cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

#kunkle='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/Kunkle_etal_Stage1_qc.gz'
kunkle='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224'
path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap'

for chr in {1..22}
do
	python aggregate_finemapper_results_modified.py \
		--out-prefix $path/max_snp_${max_snp}/kunkle \
		--sumstats $kunkle/kunkle.${chr}.snpvar_constrained.gz \
		--out $path/max_snp_${max_snp}/chr${chr}.aggregate.all.txt.gz \
		--allow-missing \
       		--chr $chr
 
done

