#!/bin/bash
for i in {1..22}
do
  echo "Running QC on chr$i"
  sbatch --export=i=$i target_job.sh
done
