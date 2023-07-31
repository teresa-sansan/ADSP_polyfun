#!/bin/bash
#SBATCH --job-name=filt_mind
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/26k_processed/%x%j.log

#21K 
#VCF_DIR="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk"
#PLINK_DIR='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/processing_files'
PLINK_DIR='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc/26k_processed'

~/plink \
  --bfile $PLINK_DIR/ADSP.geno.chr${chr}.chunk${chunk} \
  --mind 0.01 \
  --missing \
  --make-bed --out $PLINK_DIR/ADSP.mind.chr${chr}.chunk${chunk}

#done

#  --hwe 1e-6 \ ## Didn't include. Because we dont know which sample is control.
#--mind 0.01 \  ## done in vcf
#--maf 0.001 \ ## done in vcftools (set to 0.1%)


## --geno [maximum per-variant]
## --mind [maximum per-sample]

# ~/plink \
#   --bfile $PLINK_DIR/ADSP.chr${chr}.chunk${chunk} \
#   --maf 0.001 \
#   --geno 0.01 \
#   --make-bed --out $PLINK_DIR/qc/ADSP.chr${chr}.chunk${chunk}

## original bfile
## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/all/ADSP.chr${chr}.chunk${chunk}