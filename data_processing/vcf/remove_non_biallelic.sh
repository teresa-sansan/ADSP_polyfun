## a wrap script that works perfection for the 36k bigwig

for chr in {4..10}
do
    chunk_num=$(ls /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
    for ((i=1; i<=$chunk_num; i++)); 
        do
            echo check annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}
            #sbatch --mem=50G -c 10 --wrap="bcftools view /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz --regions chr${i}  > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}"
            sbatch --mem=50G -c 10 --wrap="bcftools view -m2 -M2 -v snps /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz > /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/annotated_chunk_biallelic/ADSP.chr${chr}.chunk${i}"
        done
done

