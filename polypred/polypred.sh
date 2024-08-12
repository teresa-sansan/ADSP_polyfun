#!/bin/bash
#SBATCH --job-name=bellenguez
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=70G
#SBATCH --time=5:00:00
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/%x_%j.log 

cd /gpfs/commons/home/tlin/polyfun_omer_repo
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun

##kunkle
if false; then
##susie
path='/gpfs/commons/home/tlin/output/kunkle/new_anno/'
fi

## bellenguez
if true; then
## susie
# path='/gpfs/commons/home/tlin/output/bellenguez/new_sep22/'
# max_snp=5
# anno=susie
## others
#path='/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/'
path='/gpfs/commons/home/tlin/output/prs/new_anno_0318_24/bellenguez_adsp_reference/'

#path='/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl/'
fi

# ##wightman
# if false; then
# #path='/gpfs/commons/home/tlin/output/wightman/new_anno_0203/'
# path='/gpfs/commons/home/tlin/output/wightman/new_anno_0824/'

# fi

# ## jansen
# if false; then
# path='/gpfs/commons/home/tlin/output/jansen/finemap'
# #path='/gpfs/commons/home/tlin/output/jansen/susie'
# fi

agg_file='aggregate.all.txt'

##adj beta
#agg_file='adj_beta_aggregate.all.txt'



# sconverge_file='agg_fixed_converge.tsv.gz'
python polypred.py \
	--predict \
	--betas $path/$anno/finemap/max_snp_${max_snp}/$agg_file \
	--output-prefix /gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/${anno}/finemap/polypred/36k_ibd/36k_ibd_max_snp_${max_snp}_polypred.tsv \
	--plink-exe ~/plink \
	/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg38_plink/ADSP.chr${chr}
	
#/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr*bed
	
## output /gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/${anno}/finemap/polypred/36k_ibd_max_snp_${max_snp}_polypred.tsv \
	##/gpfs/commons/home/tlin/output/wightman/new_anno_0824/${anno}/finemap/polypred/36k_ibd/36k_ibd_max_snp_${max_snp}_polypred.tsv \
	#/gpfs/commons/home/tlin/output/bellenguez/new_sep22/${anno}/finemap/polypred/36k_ibd_max_snp_${max_snp}_polypred.tsv
##	36K /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/qc/qc_chr*bed
	
	###/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/remove_triallelic/chr*_no_dups.bed
    ###17K /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/plink/vcf_filt/ADSP_annotated_merged_qc.bed
