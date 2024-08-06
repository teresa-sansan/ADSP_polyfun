#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --time=24:00:00
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err

module load snakemake

snakemake --jobs 8 --cluster "sbatch -N 1 -n {threads} --time=24:00:00"
#snakemake --jobs <number_of_jobs> --cluster "sbatch -N 1 -n <tasks_per_job> --time=<walltime>"
snakemake --jobs 5 --cluster "sbatch -N 2 -n 4 --time=3:00:00 --mem=10G"
snakemake --jobs 25 --cluster "sbatch -N 1 -n 1 --time=1:00:00 --mem=10G"