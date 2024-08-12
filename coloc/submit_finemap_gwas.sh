#!/bin/bash
#SBATCH --job-name=coloc_gwas
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --array=1-104%15


line=$SLURM_ARRAY_TASK_ID
# line=1
start=$((line +1))
head -n $start /gpfs/commons/home/tlin/script/coloc/snp_of_interest.tsv | tail -1 | while IFS=$'\t' read -r CHR anno BP
do
    name=${anno}_chr${CHR}_${BP}
    echo run $name
    sbatch --job-name=$name finemap_susie_rss.sh $CHR $BP $anno
done

