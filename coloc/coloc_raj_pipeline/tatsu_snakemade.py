import os
import numpy as np
import pandas as pd
import subprocess


dataname = "genecode_isomiga_leafcutter"
dir_base = "$working_dir"

os.chdir(dir_base)
for sub_dir in ["cluster", "tmp", "finemapping", "data"]:
	os.makedirs(os.path.join(dir_base, sub_dir), exist_ok=True)

dir_polyfun = "/sc/arion/projects/bigbrain/spliceDL/QTL/pipelines/finemapping/polyfun"
dir_geno = "/sc/arion/projects/bigbrain/spliceDL/MiGA/Genotype"


# Window for fine-mapping
window = 1000000


# Load a gene list to fine-map
# genelist = np.loadtxt("$your_path", dtype=str)

""" # If your data are sQTL and have a certain P-value threshold, a gene list can be created as below
"""

thr_p = 1e-05
top_assoc = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/microglia_eqtl_chr{0}.tsv".format(chrom), header = 0, index_col=0, sep="\t|\s+", engine="python")
genelist = top_assoc.loc[top_assoc.Random_P<thr_p].index

# Assuming a dataframe where the index is the gene name and a column 'n' represents the sample size for that gene
df_cohorts = pd.read_table("$your_path", header=0, index_col=0)



rule all:
     input:
        expand("tmp/{gene}_sumstats_hg38.gz", gene=genelist),
        expand("tmp/{gene}_sumstats.hg38.parquet", gene=genelist),
        expand("finemapping/susie.{gene}.hg38.txt.gz" % name_window, gene=genelist),
        expand("finemapping/finemap.{gene}.hg38.txt.gz" % name_window, gene=genelist),


# Extract specific gene data from full summary statistics
rule extract:
    input:
    	"data/{0}_full_assoc.tsv.gz".format(dataname)
    output:
    	sumstats_hg38="tmp/{gene}_sumstats_hg38.gz",
    shell:
    	"""
    	zcat {input} | grep 'variant_id\\|{wildcards.gene}'| gzip -f > {output.sumstats_hg38}
		"""

# If you want to use Random Z-scores, specify "--random True"
rule munge_sumstats:
    input:
	    "tmp/{gene}_sumstats_hg38.gz"
    output:
    	parquet="tmp/{gene}_sumstats.hg38.parquet",
    params:
        n=lambda wildcards: int(df_cohorts.loc[wildcards.gene, "n"]),
        #dir_polyfun=dir_polyfun
        dir_polyfun="/sc/arion/projects/bigbrain/spliceDL/conda/polyfun"
    conda:
    	"polyfun"
    shell:
	    """
	    python {params.dir_polyfun}/munge_eqtl_sumstats.revised.py --sumstats {input} --n {params.n} --out {output.parquet}
		"""

# https://github.com/omerwe/polyfun/issues/112
# If you get an error like "Estimating residual variance failed: the estimated value is negative", try specifying "--susie-resvar 1"
rule finemapping_susie:
    input:
    	sumstats="tmp/{gene}_sumstats.hg38.parquet",
    output:
    	out="finemapping/susie.{gene}.hg38.txt.gz" % name_window,
    params:
    	chr=lambda wildcards: wildcards.gene.split(":")[0].replace("chr", ""),
    	start=lambda wildcards: int(max(0, int(wildcards.gene.split(":")[1])-window)),
    	end=lambda wildcards: int(int(wildcards.gene.split(":")[2])+window),
    	dir_polyfun=dir_polyfun,
    	dir_geno=dir_geno,
    	n=lambda wildcards: int(df_cohorts.loc[wildcards.gene, "n"]),
    	geno=lambda wildcards: df_cohorts.loc[wildcards.gene, "genotype"],
    conda:
    	"polyfun3"
    shell:
    	"""
    	PYTHONPATH=/sc/arion/projects/bigbrain/spliceDL/QTL/pipelines/finemapping/polyfun
	    python {params.dir_polyfun}/finemapper.py \
	    --geno {params.dir_geno}/microglia_EUR.chr{params.chr} \
	    --sumstats {input.sumstats} --n {params.n} \
	    --chr {params.chr} --start {params.start} --end {params.end} \
	    --method susie --max-num-causal 5 \
		--out {output.out} --non-funct --allow-missing \
		--allow-swapped-indel-alleles
		#--susie-resvar 1
		"""

rule finemapping_finemap:
    input:
    	sumstats="tmp/{gene}_sumstats.hg38.parquet",
    output:
    	out="finemapping/finemap.{gene}.hg38.txt.gz" % name_window,
    params:
    	chr=lambda wildcards: wildcards.gene.split(":")[0].replace("chr", ""),
    	start=lambda wildcards: int(max(0, int(wildcards.gene.split(":")[1])-window)),
    	end=lambda wildcards: int(int(wildcards.gene.split(":")[2])+window),
    	dir_polyfun=dir_polyfun,
    	dir_geno=dir_geno,
    	n=lambda wildcards: int(df_cohorts.loc[wildcards.gene, "n"]),
    	geno=lambda wildcards: df_cohorts.loc[wildcards.gene, "genotype"],
    conda:
    	"polyfun3"
    shell:
    	"""
    	PYTHONPATH=/sc/arion/projects/bigbrain/spliceDL/QTL/pipelines/finemapping/polyfun
	    python {params.dir_polyfun}/finemapper.py \
	    --geno {params.dir_geno}/microglia_EUR.chr{params.chr} \
	    --sumstats {input.sumstats} --n {params.n} \
	    --chr {params.chr} --start {params.start} --end {params.end} \
	    --method finemap --max-num-causal 5 \
	    --finemap-exe {params.dir_polyfun}/bin/finemap_v1.4.2_x86_64 \
		--out {output.out} --non-funct --allow-missing \
		--allow-swapped-indel-alleles
		"""