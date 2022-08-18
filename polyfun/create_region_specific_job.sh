cd /gpfs/commons/home/tlin/polyfun_omer_repo
for max_snp in 1 5 10
do
	echo start creating jobs in max_snp_${max_snp} ...
	for chr in  10
	do
		python create_finemapper_jobs.py \
			--sumstats /gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl_dl_annotations/bl_dl_annotations.${chr}.snpvar_constrained.gz \
   	 		--n 63926 \
			--method susie \
    			--max-num-causal $max_snp \
			--out-prefix /gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224_annotations/bl_dl_annotations/ \
    			--jobs-file test_job.sh

	done

done
    			#--sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.${chr}.snpvar_constrained.gz \
