#!/bin/bash
#SBATCH --job-name=polyfun_pip_thres
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=15G
#SBATCH --time=10:00:00
#SBATCH --array=1-22%22
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/polyfun/pip_thres/middle_file/%x_%j.log

## calculate PRS for specific region/SNPs. So didnt perform clumping here
## q-score-range and --extract can be omitted if you are not interested in how different p value threshold perform.

polyfun_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/snp_pip_thres/'
genotype='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.chr'
output='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/polyfun/pip_thres/middle_file/'
sumstat='bellenguez'
chr=$SLURM_ARRAY_TASK_ID

for anno in baseline omics omics_dl
do 
    for thres in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 
    do

    if [ "$anno" != "susie" ]; then 
        echo running $anno
            ~/plink \
            --bfile ${genotype}${chr} \
            --extract $polyfun_path/${sumstat}_${anno}_pip_$thres.snp \
            --score $polyfun_path/remove_index0/${sumstat}_${anno}.txt 2 4 7 header \
            --out $output/${anno}_thres${thres}_chr${chr}
    else
        echo running $anno
            ~/plink \
            --bfile ${genotype}${chr} \
            --extract $polyfun_path/${sumstat}_${anno}_pip_$thres.snp \
            --score $polyfun_path/remove_index0/${sumstat}_${anno}.txt 3 5 6 header \
            --out $output/${anno}_thres${thres}_chr${chr}
        fi
    done
done

#/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/${sumstats}.snp 
#baseline omics omics_dl 

## susie : 1,4,11 ## also snp and chr is reversed! 
## others: 2,4,12