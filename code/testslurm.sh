#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --job-name=testjob
#SBATCH --output=testslurm.{file}.out
echo "hello"

echo ${file}


