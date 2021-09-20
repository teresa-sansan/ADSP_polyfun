#!/bin/bash
outputformat='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/bellenguez_all'
#outputformat='/gpfs/commons/home/tlin/output/kunkle_all_2/all_anno'
#outputformat=/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_deepsea_brain_atac/bl_roadmap_deepsea_brain_atac
#outputformat=/gpfs/commons/home/tlin/polyfun/output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac
#outputformat=output/bl_roadmap_deepsea/bl_roadmap_deepsea
#outputformat=output/jansenetal/bl_jansen
#outputformat=output/kunkle_baseline_microglia_atac/kunkle_bl_microglia
#outputformat=output/bl_roadmap/bl_roadmap
#outputformat=output/bl_roadmap/specific_col/bl_roadmap
#outputformat=output/bl_roadmap_microglia/bl_roadmap_microglia
#outputformat=output/bl_deepsea/specific_col/bl_deepsea


for i in {12..22}
do

sbatch --export=chr=$i,output=$outputformat polyfun_1_3.sh

done


#sbatch --export=chr=$i,output=$outputformat/H3K27me3/bl_deepsea_H3K27me3_ polyfun_1_3.sh
#sbatch --export=chr=$i,output=$outputformat/H3K4me3/bl_deepsea_H3K4me3_ polyfun_1_3.sh
#sbatch --export=chr=$i,output=$outputformat/H3K4me1/bl_deepsea_H3K4me1_ polyfun_1_3.sh
#sbatch --export=chr=$i,output=$outputformat/H3K9ac/bl_deepsea_H3K9ac_ polyfun_1_3.sh
#sbatch --export=chr=$i,output=$outputformat/H3K27ac/bl_deepsea_H3K27ac_ polyfun_1_3.sh
