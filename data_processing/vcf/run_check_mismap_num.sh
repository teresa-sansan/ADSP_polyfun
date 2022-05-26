#!/bin/bash

for chr in {1..22}
do
	sbatch --export=i=$chr check_mismap_num.sh
done
