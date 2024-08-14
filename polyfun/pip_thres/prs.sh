#!/bin/bash
#SBATCH --job-name=prs_polyfun
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=16G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa/prs/%x_%j.log

dir_beta="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa/"
dir_geno="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc"

# anno=$1
# if [ $anno == "susie" ]; then
#     snp=3
#     A1=5
#     beta=11
# else
#     snp=2
#     A1=4
#     beta=12
# fi

# for thres in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
# do
#     for chr in {1..22}
#     do
#         ~/plink \
#             --bfile  $dir_geno/ADSP.chr${chr} \
#             --score $dir_beta/bellenguez_remove_index0_${anno}.txt $snp $A1 $beta header \
#             --extract $dir_beta/pip_thres/${anno}_pip${thres}.snp \
#             --out $dir_beta/prs/middle_file/${anno}_pip${thres}.chr${chr}
#     done
#     awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $dir_beta/prs/middle_file/${anno}_pip${thres}.chr*.profile > $dir_beta/prs/middle_file/${anno}_pip_${thres}.prs

# done 

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun
python merge_diagnosis_toplot.py $anno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa/prs/

