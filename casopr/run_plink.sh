#!/bin/bash
#SBATCH --job-name=casioPR_prs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/pic/casioPR/simulation/prs/0214_no_vld_inf/%x%j.log

path='/gpfs/commons/home/tlin/pic/casioPR/simulation/'
wightman_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp'


#!/bin/bash


#Function to run PLINK command
run_plink() {
    local chr="$1"
    local score="$2"
    local iter=$((score - 12))

    path='/gpfs/commons/home/tlin/pic/casioPR/simulation/'
    #path='/gpfs/commons/home/tlin/pic/casioPR/simulation/prs/intersect_2'
    wightman_snp='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/wightman_fixed_beta.snp'
    plink_path="/gpfs/commons/home/tlin/plink"
    #prs_path="$path/prs/0128_chr1_22_iter5_v0e_betas.tsv"
    prs_path="$path/prs/intersect_2/0219_no_vld_inf_chr1_22_iter3_jab_betas.tsv"
    #prs_path='/gpfs/commons/home/tlin/pic/casioPR/simulation/prs/0214_no_vld_inf/0214_no_vld_inf_chr1_22_iter10_mmv_betas.tsv'
    name='0219_inf_whole_genome_jab_iter'
    range_file="/gpfs/commons/home/tlin/script/plink/range_list.txt"

    /gpfs/commons/home/tlin/plink \
    --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr${chr} \
    --score $prs_path 2 4 $score \
    --q-score-range $range_file $wightman_snp \
    --out $path/prs/intersect_2/${name}iter${iter}_prs_chr${chr}


    #--out $path/prs/0214_no_vld_inf/0214_inf_whole_genome_mmv_iter${iter}_prs_chr${chr}
    
}


# Run PLINK commands in parallel for chr 1 to 22 and for each score value
export -f run_plink 
parallel -j 3 run_plink ::: {13..22} ::: {12..14} 


# for chr in 1
# do
# ~/plink \
# --bfile /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr$chr \
# --score /gpfs/commons/home/tlin/pic/casioPR/simulation/prs/0128_chr1_22_iter5_v0e_betas.tsv 2 4 13 \
# --q-score-range /gpfs/commons/home/tlin/script/plink/range_list2.txt $wightman_snp \
# --out $path/prs/bl_0128_chr1_22_iter5_v0e_iter1_prs_chr$chr
# done


#score: need to followed by SNP ID, effective allele, effect size


