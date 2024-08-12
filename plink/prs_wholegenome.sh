#!/bin/bash
#SBATCH --job-name=bellenguez_adsp_reference
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=7:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/v2/%x_%j.log

## calculate PRS for specific region/SNPs. So didnt perform clumping here
## Note: only need to calculate the chr that the SNPs of interests located. ## (i.e. for APOE, you only need to look into chr19)
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

polyfun_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap_v2/combined'
output='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/v2'

#chr=$SLURM_ARRAY_TASK_ID
anno=$1
echo running $anno
if [ "$anno" == "susie" ]; then 
  SNP=1
  beta_mean=11
else 
  SNP=2
  beta_mean=12
fi

zcat $polyfun_path/bellenguez_${anno}.txt.gz > $output/bellenguez_${anno}.txt
for chr in {1..22}
do
    ~/plink \
    --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink/ADSP.chr${chr} \
    --score  $output/bellenguez_${anno}.txt $SNP 5 $beta_mean header \
    --out $output/${anno}_chr${chr}

done
rm $output/bellenguez_adsp_reference/${anno}_chr${chr}.txt


awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $output/${anno}_chr*.profile > $output/${anno}.prs
echo wrote $output/${anno}.prs

