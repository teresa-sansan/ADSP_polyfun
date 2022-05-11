cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37

zcat ADSP_annotated_chr3.vcf.gz| awk 'NF == 7 {print; exit 1}'
