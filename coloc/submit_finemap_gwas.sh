#!/bin/bash
#SBATCH --job-name=coloc_gwas
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=30:00:00
#SBATCH --array=1-324%50

#/gpfs/commons/home/tlin/script/coloc/snp_of_interest.tsv
line=$SLURM_ARRAY_TASK_ID
# line=1
start=$((line +1))
head -n $start /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/snp_pip_thres/remove_index0/union_snp_of_interest.tsv | tail -1 | while IFS=$'\t' read -r CHR anno BP
do
    name=${anno}_chr${CHR}_${BP}
    echo run $name
    sbatch --job-name=$name finemap_susie_rss.sh $CHR $BP $anno
done

