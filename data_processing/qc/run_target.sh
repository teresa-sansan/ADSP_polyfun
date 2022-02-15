#!/bin/bash
for i in {1..19}
do
  echo "Running QC on chr$i"
  sbatch --export=i=$i target_job.sh
done
