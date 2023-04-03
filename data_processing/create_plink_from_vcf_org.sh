path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37/'

for i in {1..14}   ##1..14
do
echo running chr $i
sbatch --mem=90G -c 5 --wrap="~/plink --vcf $path/correct_chr_vcf/filt/ADSP_annotated_chr${i}.genofilt.vcf.recode.vcf --make-bed --double-id -allow-extra-chr --out $path/plink/vcf_filt/ADSP_annotated_chr$i"
done


## added two flags (--double-id & --allow_extra-chr) because
## there were 
## Error: Multiple instances of '_' in sample ID. (If you do not want '_' to be treated as a FID/IID delimiter, use --double-id )
## Error: Invalid chromosome code 'chr17_gl000205_random' on line 1790 of .vcf file. (e.g. in chr21)
## increase the memrory size. ( bus error)


##zero-cluster flag doesnt work
##--zero-cluster plink/zero_cluster_chr$i
##Error: --zero-cluster must be used with --within/--family.
