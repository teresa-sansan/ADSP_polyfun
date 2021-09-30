cd /gpfs/commons/home/tlin/polyfun_omer_repo
chr=22
python polypred.py \
	--predict \
	--betas /gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/chr${chr}.aggregrate.all.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/kunkle_all/polypred/polypred_kunkle_all2_chr${chr} \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/filtered_data/ADSP_chr*.bed
