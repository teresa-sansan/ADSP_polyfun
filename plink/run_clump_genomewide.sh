#!/bin/bash

## no qc
#sbatch --export=sumstats=kunkle,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Kunkle_et_al_2019_hg37_ldsc.tsv.gz' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh
#sbatch --export=sumstats=bellenguez,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Bellenguez_et_al_2021_hg37_ldsc.tsv.gz' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh
#sbatch --export=sumstats=wightman,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Wightman_et_al_2021_hg37_ldsc.tsv.gz' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh


# qc
#sbatch --export=sumstats=kunkle,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Kunkle_et_al_2019_hg37_ldsc_qc.tsv' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh
#sbatch --export=sumstats=bellenguez,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_qc.tsv' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh
sbatch --export=sumstats=wightman,sumstats_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_qc.tsv' /gpfs/commons/home/tlin/script/plink/clump_genomewide.sh
