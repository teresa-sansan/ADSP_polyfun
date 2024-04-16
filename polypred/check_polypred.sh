python /gpfs/commons/home/tlin/polyfun_omer_repo/polypred.py \
	--predict \
	--betas /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/bellenguez_baseline_chr19.txt.gz \
	--output-prefix /gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/check_result_bl_19.tsv \
	--plink-exe /gpfs/commons/home/tlin/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink/ADSP.chr19.bed
