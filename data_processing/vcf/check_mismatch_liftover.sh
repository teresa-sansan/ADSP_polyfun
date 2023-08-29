#!/bin/bash
#SBATCH --job-name=liftover_check
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G 
#SBATCH --time=100:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37/%x%j.log

i=22

cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37/
#zcat ADSP_chr${i}.vcf.gz| grep -v "##"| cut -f 1| uniq > chr${i}_uniq2.txt
zcat ADSP_chr${i}.vcf.gz| grep -v "##"| cut -f 1| sort | uniq -c > chr${i}_check2.txt

# cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated
# touch org_check_2.txt
# for i in {1..22}
# do
# echo chr${i}
# zcat ADSP.chr${i}.vcf.gz|grep -v "##"| cut -f 1| wc -l >> org_check2.txt 
# done