#!/bin/bash
#SBATCH --job-name=prs_bellenguez_remove_index0
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --array=1-22%3
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/remove_index0/%x_%j.log

## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

polyfun_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fine_mapping/annotations_dl/aggregate_finemap/'
adsp_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
#output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/remove_index0/'
sumstat='bellenguez'
thres=$1
chr=$SLURM_ARRAY_TASK_ID

for anno in baseline omics omics_dl susie 
do 
  echo running $anno
  #zcat $adsp_path/${sumstat}_${anno}_chr${chr}.txt.gz > $output/bellenguez_adsp_reference/${sumstat}_${anno}_chr${chr}.txt
  if [ "$anno" == "susie" ]; then 
    SNP=1
    beta_mean=11
  else 
    SNP=2
    beta_mean=12
  fi

  ~/plink \
  --allow-no-sex \
  --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink/ADSP.chr${chr} \
  --extract $adsp_path/pip_filt/${sumstat}_${anno}_pip_${thres}.snp \
  --score  $adsp_path/remove_index0/${sumstat}_${anno}_chr${chr}.txt $SNP 5 $beta_mean header \
  --out $output/${anno}_chr${chr}_pip_${thres} 

  # rm  $output/bellenguez_adsp_reference/allele_flip/${anno}_chr${chr}.nosex
  # rm  $output/bellenguez_adsp_reference/${sumstat}_${anno}_chr${chr}.txt
done

#baseline omics omics_dl  susie
#  --extract $adsp_path/pip_filt/${sumstat}_${anno}_pip_${thres}.snp \
 # --extract /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/pip_filt/${sumstat}_${anno}_pip_${thres}.snp \
## susie : 1,4,11 ## also snp and chr is reversed! 
## others: 2,4,12
## old /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr
## --out $output/bellenguez_adsp_reference/thres/plink_output/sum_${anno}_pip_${thres}_chr${chr} 

#--score  $output/bellenguez_adsp_reference/${sumstat}_${anno}_chr${chr}.txt $SNP 5 $beta_mean header \
#