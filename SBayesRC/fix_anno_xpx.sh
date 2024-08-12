path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC'
# chrom=22
# echo $path/fix_chr$chrom.tsv 
# #cat $path/chr${chrom}.tsv | cut -f 109,110 --complement > $path/fix_chr${chrom}.tsv
# cat $path/fix_chr${chrom}.tsv | cut -f 92,91,90 --complement > $path/fixed_chr${chrom}.tsv
# mv $path/fixed_chr${chrom}.tsv $path/fix_chr${chrom}.tsv
# head -1 $path/fix_chr${chrom}.tsv | wc



#cat $path/new_all.tsv | cut -f 110 --complement > $path/fix_all.tsv
cat $path/fix_all.tsv | cut -f 96 --complement > $path/fixed_all.tsv
mv $path/fixed_all.tsv $path/fix_all.tsv
head -1 $path/fix_all.tsv | wc
