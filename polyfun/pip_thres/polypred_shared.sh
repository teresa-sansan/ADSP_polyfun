#!/bin/bash
#SBATCH --job-name=bellenguez
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=30G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

## ADSP_36K_filtered
pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')

## bellenguez
path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/'
agg_file='aggregate.all.txt'

python polypred.py \
	--predict \
	--betas $path/$anno/finemap/max_snp_${max_snp}/$agg_file \
	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/${anno}/finemap/polypred/36k_ibd/36k_ibd_max_snp_${max_snp}.prs \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink_qc/ADSP.*bed
	


