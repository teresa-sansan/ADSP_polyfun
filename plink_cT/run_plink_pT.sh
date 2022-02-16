#!/bin/bash
#echo creat snp_pvalue file
#zcat /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/bellenguez_2021_final_rename.tsv.gz| awk '{print $3,$9}' > /gpfs/commons/home/tlin/data/bellenguez_snp.pvalue

for i in {1..22}
do
  echo "start calculating PRS for chr $i"
  sbatch --export=chr=$i kunkle_plink_pT.sh   
  #sbatch --export=chr=$i /gpfs/commons/home/tlin/polyfun_script/test/pT_height.sh
done
