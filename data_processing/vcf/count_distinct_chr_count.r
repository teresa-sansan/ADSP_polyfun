library(dplyr)

for (i in 1:22) {
print(paste("starting chr",i ))
file_name=(paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/incorrect_vcf/ADSP_annotated_hg37_chr',toString(i),'_incorrect.vcf.tsv', sep = ''))
output_file=(paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/only_position/incorrect_vcf/wrong_mapping_chr_count_chr',toString(i),'.tsv', sep = ''))

file=read.csv(file_name, sep = '\t')
count_wrong_chr = file %>%
	count(file$X.CHROM)

write.table(count_wrong_chr, output_file, row.names=FALSE, col.names=c("chr","count"),sep = '\t')
}
