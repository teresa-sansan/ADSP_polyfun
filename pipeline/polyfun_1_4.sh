#!/bin/bash
#SBATCH --job-name=polyfun1-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=output/bellenguez/bellenguez_bl/%x_%j.log

cd ~/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
output='/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/bellenguez_bl'
sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'

python polyfun.py \
 	--compute-h2-bins \
    	--output-prefix $output \
    	--sumstats $sumstat \
    	--w-ld-chr $bl/weights.UKB. \
    	--allow-missing

#output='/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_deepsea_brain_atac/bl_roadmap_deepsea_brain_atac'
#output='/gpfs/commons/home/tlin/polyfun/output/bl_microglia/kunkle_bl_microglia'
#output='/gpfs/commons/home/tlin/polyfun/output/bl_deepsea_microglia/bl_deepsea_microglia' 
#output='/gpfs/commons/home/tlin/polyfun/output/bl_deepsea_microglia/bl_deepsea_microglia'
#output='/gpfs/commons/home/tlin/polyfun/output/bellenguez/bellenguez/bellenguez'
#output='output/bl_roadmap_microglia/bl_roadmap_microglia'
#output='/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_microglia/bl_roadmap_microglia'
#sumstat='/gfps/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/jansenetal_2019/AD_sumstats_Jansenetal_2019sept.parquet'
#output=output/jansenetal/bl_jansen
#output='output/bl_roadmap/all_col/bl_roadmap'
#sumstat='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet'
#output='/gpfs/commons/home/tlin/polyfun/output/bl_deepsea/specific_col/bl_deepsea'

#for i in H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9ac
#do
#	python polyfun.py \	
#    		--compute-h2-bins \
#    		--output-prefix /gpfs/commons/home/tlin/polyfun/output/bl_deepsea/$i/bl_deepsea_${i}_ \
#    		--sumstats $sumstat \
#    		--w-ld-chr $bl/weights.UKB. \
#    		--allow-missing
#done
