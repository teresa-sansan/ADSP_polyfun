#!/bin/bash
#SBATCH --job-name=run_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=120G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/home/tlin/polyfun_script/pipeline/create_job/%x%j.log


for chr in {1..21}
do
	bash finemap_all_jobs_secondary.${chr}.txt 
	echo submit chr $chr
done
