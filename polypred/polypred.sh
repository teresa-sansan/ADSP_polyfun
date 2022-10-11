#!/bin/bash
#SBATCH --job-name=Polypred_bellenguez_not_fixed
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap/polypred_new_plink/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

##kunkle
if false; then
##susie
path='/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/susie_finemap'
fi

## bellenguez
if true; then
## susie
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/susie'
path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224/finemap'
#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl/'
fi

##wightman
if false; then
#path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/all_anno/finemap/'
path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/'
fi

## jansen
if false; then
path='/gpfs/commons/home/tlin/output/jansen/finemap'
#path='/gpfs/commons/home/tlin/output/jansen/susie'
fi

python polypred.py \
	--predict \
	--betas $path/max_snp_${max_snp}/aggregate.all.txt \
	--output-prefix $path/polypred_new_plink/not_fixed_max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_*.bed


