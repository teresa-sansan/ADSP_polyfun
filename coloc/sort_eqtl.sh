#!/bin/bash
#SBATCH --job-name=eqtl
#SBATCH --partition=pe2
#SBATCH --nodes=1           # minimum number of nodes to be allocated
#SBATCH --ntasks=1          # number of tasks
#SBATCH --cpus-per-task=8   # number of cores on the CPU for the task
#SBATCH --mem=15G
#SBATCH --time=3:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --array=18-22%5
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/%x_%j.log

chr=chr$SLURM_ARRAY_TASK_ID
echo $chr
echo 'sorting..'

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia

cat header.txt > eqtl_$chr.tsv
zcat ../microglia_GENCODE_expression_full_assoc.tsv.gz | awk -v chr=$chr '{ if($3 == chr) {print}}' > eqtl_$chr.tsv

# done
# head -n 3 eqtl_$chr.tsv
# wc -l eqtl_$chr.tsv