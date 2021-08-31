cd /gpfs/commons/home/tlin/polyfun_omer_repo

for chr in {1..22}
do
	python create_finemapper_jobs_modified.py \
    		--sumstats /gpfs/commons/home/tlin/output/kunkle_all/all_anno_.${chr}.snpvar_constrained.gz \
   	 	--n 63926 \
		--method susie \
    		--max-num-causal 1 \
    		--out-prefix /gpfs/commons/home/tlin/output/kunkle_all/finemap_genomewide/all_anno \
    		--jobs-file /gpfs/commons/home/tlin/polyfun_script/pipeline/create_job/finemap_all_jobs.${chr}.txt

done
