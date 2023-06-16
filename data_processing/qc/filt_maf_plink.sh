#!/bin/bash
#SBATCH --job-name=filt_mind
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=22:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38/qc/%x%j.log


VCF_DIR="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk"
PLINK_DIR='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38'

~/plink \
  --bfile $PLINK_DIR/qc/ADSP.chr${chr}.chunk${chunk} \
  --mind 0.01 \
  --make-bed --out $PLINK_DIR/qc/ADSP.mind.chr${chr}.chunk${chunk}

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