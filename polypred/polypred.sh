#!/bin/bash
#SBATCH --job-name=Polypred
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=70G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/new_anno_0203/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

##kunkle
if false; then
##susie
path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap'
fi

## bellenguez
if false; then
## susie
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/susie'
path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl/'
fi

##wightman
if true; then
path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/'
fi

## jansen
if false; then
path='/gpfs/commons/home/tlin/output/jansen/finemap'
#path='/gpfs/commons/home/tlin/output/jansen/susie'
fi

agg_file='aggregate.all.txt'
# sconverge_file='agg_fixed_converge.tsv.gz'

python polypred.py \
	--predict \
	--betas $path/$anno/finemap/max_snp_${max_snp}/$agg_file \
	--output-prefix $path/$anno/finemap/polypred/max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc.bed
	
	
	###/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/remove_triallelic/chr*_no_dups.bed

