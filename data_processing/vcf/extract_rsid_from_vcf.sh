#!/bin/bash
#SBATCH --job-name=36k_RSID
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=100:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/PRS_hg38/vcf_rsid/%x_%j.log


cd /gpfs/commons/home/tlin/PRScs
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_preview/'
# chunk_num=$(ls $path/plink_hg38| grep chr"${chr}".chunk | grep bed | wc -l)
# echo chumk num is $chunk_num

#for ((i=1; i<=$chunk_num; i++)); 
#do
#zcat $path/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz|grep -v "#"|cut -f 3| grep -v '.' > $path/PRS_hg38/vcf_rsid/ADSP.chr${chr}.chunk${i}.snp
#tabix -h -R $region $path/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz | grep -v "#" | cut -f 3 | grep -v '.' > $path/PRS_hg38/vcf_rsid/ADSP.chr${chr}.chunk${i}.snp
tabix -h $path/annotated_chunk/ADSP.chr${chr}.chunk${i}.vcf.bgz chr${chr} | grep -v "#" | cut -f 3 | grep -v '.' > $path/PRS_hg38/vcf_rsid/new_ADSP.chr${chr}.chunk${i}.snp


#done

