#!/bin/bash
#SBATCH --job-name=parallel_36k_ibd
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=35G
#SBATCH --time=60:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez/%x_%j.log

#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/pre_qc'

# for chr in {21..22}
# do
# echo starting chr$chr ...
# cat $path/ADSP.geno.chr${chr}.bim| cut -f 2 | sort | uniq -d > ${chr}.dup
# done

# Define the function to process each chromosome
process_chr() {
    path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/pre_qc'
    local chr="$1"
    #cat "$path/ADSP.geno.chr${chr}.bim" | cut -f 2 | sort | uniq -d > "$path/chr${chr}.dup"
    ~/plink --bfile $path/ADSP.geno.chr${chr} --exclude $path/chr${chr}.dup --make-bed --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/pa_qc_chr${chr}

}
# Export the process_chr function so it's visible to parallel
#export -f process_chr


clump_chr(){
    plink_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc'
    bfile='qc_chr'
    local chr="$1"
    ~/plink \
    --bfile $plink_path/$bfile${chr}\
    --clump-p1 1 \
    --clump-r2 0.1  \
    --clump-kb 250  \
    --clump /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/Bellenguez_et_al_2021_hg37_new_sep20_qc.tsv \
    --clump-snp-field SNP \
    --clump-field P \
    --out /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_bellenguez/ADSP_qc_plink_${chr} 

}

export -f clump_chr

# Generate a sequence of chromosome numbers from 1 to 22
chromosomes=$(seq 9 12)

# Use parallel to process the chromosomes in parallel
#echo "$chromosomes" | parallel -j $(nproc) --eta process_chr {} --bibtex

#!/bin/bash
# Use parallel to run your script for each chromosome in parallel
echo "$chromosomes" | parallel -j $(nproc) --eta clump_chr {}
