#!/bin/bash
#SBATCH --job-name=hg38_plink
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=15:00:00
#SBATCH --array=1-22%15
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/%x_%j.log

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun 
## keep = SUBJID from /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_redo_his.tsv
## note the maf is 0.1%

chr=$SLURM_ARRAY_TASK_ID

# create hg38 plink (for bellenguez in ADSPpanel)
#variant_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/ADSP_EUR/plink/'
variant_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/ADSP_EUR/' ## this is updated in June
#output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/plink_file/' ## this is in may
output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/plink_file_hg38/' 
echo start chr $chr ...

# cd $output_path
# cut -f 2 ADSP.chr${chr}.bim | sort | uniq -d > ${chr}_dups
# ~/plink --bfile ADSP.chr${chr} --exclude ${chr}_dups --make-bed --out ADSP_chr${chr}_filt --snps-only 


#--maf 0.001 \


#plink --vcf "$1" --out "$OUT$2" --const-fid --maf 0.000000000000001 --make-bed --keep "$EXTRACT"keep.txt
#--allow-extra-chr
#--keep $variant_path/ADSP_bellenguez_snps_chr${chr}.txt \


#awk '{print $2}' $variant_path/ADSP_bellenguez_snps_chr${chr}.bim  > $variant_path/ADSP_bellenguez_snps_chr${chr}.txt

~/plink2/plink2 \
--const-fid \
--make-bed --snps-only --exclude exclude_snps.txt \
--out $output_path/ADSP_EUR_chr${chr} \
--vcf $variant_path/ADSP_EUR_chr${chr}.vcf.gz

#rm $variant_path/ADSP_bellenguez_snps_chr${chr}.txt
