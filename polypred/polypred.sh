#!/bin/bash
#SBATCH --job-name=Polypred_susie
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/polypred/%x_%j.log 

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
path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/susie/finemap/'
#path='/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/bl/finemap/'

fi

## jansen
if false; then
path='/gpfs/commons/home/tlin/output/jansen/finemap'
#path='/gpfs/commons/home/tlin/output/jansen/susie'
fi


ori_file='aggregate.all.txt'
converge_file='agg_fixed_converge.tsv.gz'

python polypred.py \
	--predict \
	--betas $path/max_snp_${max_snp}/$converge_file \
	--output-prefix $path/polypred/fixed_max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/ADSP_qc_all/ADSP_qc_all_*.bed


