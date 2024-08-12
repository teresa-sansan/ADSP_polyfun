python /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_susie_rss.py \
	--geno /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/plink_ADSP_filtered_deduplicated/ADSP_chr8 \
	--sumstats /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_backup_teresa/bellenguez_omics/bellenguez_omics.8.snpvar_ridge.gz \
	--n 487511 --chr 8 --start 143500001 --end 145500001 \
	--method susie --max-num-causal 10 \
	--no-sort-pip --susie-resvar 1 \
	--susie-outfile /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/snp_of_interest/omics_chr8_144500001 \
	--allow-missing --out /gpfs/commons/home/tlin/output/36k/bellenguez/adsp_ld/susie_rss/snp_of_interest/omics_chr8_144500001
