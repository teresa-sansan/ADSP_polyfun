#!/bin/bash
#SBATCH --job-name=run_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=120G
#SBATCH --time=19:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/kunkle_all_2/finemap_max_snp_1/%x%j.log


#for chr in {1..21}
#do
#	bash finemap_all_jobs_secondary.${chr}.txt 
#	echo submit chr $chr
#done

sbatch finemap_all_jobs_kunkle2.1.txt
