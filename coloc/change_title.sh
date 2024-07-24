cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia

chr=12
## change title
sed -e '1s/chr/CHR/' -e '1s/pos/BP/'  -e '1s/variant_id/SNP/' -e '1s/ref/A1/'  -e '1s/alt/A2/'  -e '1s/z_tissue_2/Z/' chr12.tsv > chr12_z2.tsv