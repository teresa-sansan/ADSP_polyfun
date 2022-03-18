cd /gpfs/commons/home/tlin/polyfun_omer_repo
for max_snp in  5
do
	echo start creating jobs in max_snp_${max_snp} ...
	for chr in  1	
	do
		python create_finemapper_jobs.py \
			--sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_test/bellenguez.${chr}.snpvar_constrained.gz \
   	 		--n 487511 \
			--method susie \
    			--max-num-causal $max_snp \
			--out-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224_test/finemap/max_snp_${max_snp}/bellenguez \
			--jobs-file /gpfs/commons/home/tlin/script/polyfun/finemap_job/bellenguez/check/bellenguez_max_snp_${max_snp}_chr${chr}.sh 
    			#--out-prefix /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap_fixed_assertion_susie_iter/max_snp_${max_snp}/bellenguez \
    			#--jobs-file /gpfs/commons/home/tlin/script/polyfun/finemap_job/bellenguez_max_snp_${max_snp}_chr${chr}.sh &

	done

done
    			#--sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.${chr}.snpvar_constrained.gz \
