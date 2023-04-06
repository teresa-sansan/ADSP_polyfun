#!/bin/bash 
for i in {1..19}
do
echo run chr $i
sbatch --export=chr=$i recover_header.sh
#remove_snp_vcf.sh
# compress.sh
done