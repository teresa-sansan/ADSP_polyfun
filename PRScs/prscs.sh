#!/bin/bash
#SBATCH --job-name=bellenguez_prscs_17k
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/prs/PRSCS/17k/bellenguez/%x_%j.log


cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

sumstat_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed'
bellenguez='bellenguez_4prscs.tsv'
n_bellenguez=487511
##bellenguez


for chr in {11.22}
do
python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
    --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
    --sst_file=$sumstat_path/$bellenguez \
    --n_gwas=$n_bellenguez --phi=1e-2 \
    --chrom=$chr \
    --out_dir=/gpfs/commons/home/tlin/output/prs/PRSCS/17k/bellenguez/
done

# python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
#     --bim_prefix=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc \
#     --sst_file=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_4prscs_not0.tsv \
#     --n_gwas=762971 --phi=1e-2 \
#     --out_dir=/gpfs/commons/home/tlin/output/wightman/prscs/all_anno/PIP_not0/ 



path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38_qc'
# chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
# echo chumk num is $chunk_num


# for chr in {1..22}
# do
# echo $path/ADSP.chr${chr}.chunk${i} 

# python PRScs.py --ref_dir=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_PRScs/ldblk_ukbb_eur \
#    --bim_prefix=$path/ADSP_qc_chr${chr} \
#    --sst_file=$sumstat_path/$bellenguez \
#    --chrom=$chr --n_gwas=$n_bellenguez --phi=1e-2 --out_dir=/gpfs/commons/home/tlin/output/prs/PRSCS/36k/bellenguez/chr${chr}
# done

# cat $path/PRS/prscs/chr${chr}_chunk* > $path/PRS/prscs/chr${chr}.beta 
# echo "concat all chunks in chr " $chr






## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_prscs.tsv
## --ref_dir: LD panel file
## bim_prefix: bim file for target data
## sst_file: sumstat_file 

## /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc