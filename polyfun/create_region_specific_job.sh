cd /gpfs/commons/home/tlin/polyfun_omer_repo

for chr in 1
do
	python create_finemapper_jobs_modified.py \
    		--sumstats /gpfs/commons/home/tlin/output/kunkle_all_2/all_anno.${chr}.snpvar_constrained.gz \
   	 	--n 63926 \
		--method susie \
    		--max-num-causal 1 \
    		--out-prefix /gpfs/commons/home/tlin/output/kunkle_all_2/finemap/all_anno_2 \
    		--jobs-file /gpfs/commons/home/tlin/polyfun_script/pipeline/create_job/finemap_all_jobs_kunkle2.${chr}.txt

done
