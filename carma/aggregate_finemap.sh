#!/bin/sh
#SBATCH --job-name=carma_agg
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=60:00:00
#SBATCH --array=1-21%21
#SBATCH --output=/gpfs/commons/home/tlin/output/CARMA/geno_filt/%x_%j.log

# chr=$SLURM_ARRAY_TASK_ID
chr=22
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt
n_ld=$(ls -1 | grep .bim| grep chr${chr}_ | wc -l)
echo chr $chr, total ld = $n_ld

cd /gpfs/commons/home/tlin/output/CARMA/geno_filt 
for ld in $(seq 1 "$n_ld")
do  
    echo $ld
    if [ $ld = 1 ]; then
        zcat ${chr}_1.txt.gz|awk '($12 > 0){print $0}' > finemap_chr${chr}.txt
    else
        zcat ${chr}_${ld}.txt.gz|tail -n+2| awk '($12 > 0){print $0}' >> finemap_chr${chr}.txt
    fi
    
done
