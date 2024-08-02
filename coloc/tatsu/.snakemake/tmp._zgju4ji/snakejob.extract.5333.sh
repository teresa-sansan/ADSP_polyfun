#!/bin/sh
# properties = {"type": "single", "rule": "extract", "local": false, "input": ["/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia_GENCODE_expression_full_assoc.tsv.gz"], "output": ["/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene/ENSG00000139537.11_sumstats_hg38.gz"], "wildcards": {"gene": "ENSG00000139537.11"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/scratch"}, "jobid": 5333, "cluster": {}}
 cd /gpfs/commons/home/tlin/script/coloc_raj_pipeline/tatsu && \
/gpfs/commons/home/tlin/.conda/envs/polyfun/bin/python3.6 \
-m snakemake /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene/ENSG00000139537.11_sumstats_hg38.gz --snakefile /gpfs/commons/home/tlin/script/coloc_raj_pipeline/tatsu/Snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/gpfs/commons/home/tlin/script/coloc_raj_pipeline/tatsu/.snakemake/tmp._zgju4ji' '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia_GENCODE_expression_full_assoc.tsv.gz' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules extract --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /gpfs/commons/home/tlin/.conda/envs/polyfun/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /gpfs/commons/home/tlin/script/coloc_raj_pipeline/tatsu/.snakemake/tmp._zgju4ji/5333.jobfinished || (touch /gpfs/commons/home/tlin/script/coloc_raj_pipeline/tatsu/.snakemake/tmp._zgju4ji/5333.jobfailed; exit 1)

