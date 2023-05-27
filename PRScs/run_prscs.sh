#!/bin/bash
for chr in {20..22} 
do
	sbatch --export=chr=$chr prscs.sh
done