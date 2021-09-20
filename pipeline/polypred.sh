cd /gpfs/commons/home/tlin/polyfun_omer_repo

python polypred.py \
	--predict \
	--betas /gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_1/aggregrate.all.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/kunkle_all/polypred/polypored_kunkle_all \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/filtered_data/ADSP_chr*.bed
