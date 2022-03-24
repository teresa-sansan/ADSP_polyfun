cd /gpfs/commons/home/tlin/polyfun_omer_repo
for max_snp in 1 3 5 7
do
	echo start creating jobs in max_snp_${max_snp} ...
	for chr in {1..22}
	do
		python create_finemapper_jobs.py \
    			--sumstats /gpfs/commons/home/tlin/output/wightman/wightman_all.${chr}.snpvar_constrained.gz \
   	 		--n 74004 \
			--method susie \
    			--max-num-causal $max_snp \
    			--out-prefix /gpfs/commons/home/tlin/output/wightman/finemap_fixed_assertion_susie_iter/max_snp_${max_snp}/finemap_wightman \
    			--jobs-file /gpfs/commons/home/tlin/script/polyfun/finemap_job/wightman/wightman_max_snp_${max_snp}_chr${chr}.sh &

	done

done
