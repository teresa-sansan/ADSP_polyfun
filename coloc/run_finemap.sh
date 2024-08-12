#sbatch --job-name="finemap_susie_chr12" ../polyfun/finemap_susie_rss.sh 12 37371914 susie
#sbatch --job-name="finemap_susie_chr8" ../polyfun/finemap_susie_rss.sh 8 142943835 susie
#sbatch --job-name="finemap_susie_chr21_susie" ../polyfun/finemap_susie_rss.sh 21 34567305 susie
# sbatch --job-name="finemap_susie_chr21_ENSG00000159140.21" ../polyfun/finemap_susie_rss.sh 21 34567305 ENSG00000159140.21
# sbatch --job-name="finemap_susie_chr21_ENSG00000234380.2" ../polyfun/finemap_susie_rss.sh 21 34567305 ENSG00000234380.2
# sbatch --job-name="finemap_susie_chr21_ENSG00000249624.10" ../polyfun/finemap_susie_rss.sh 21 34567305 ENSG00000249624.10
# sbatch --job-name="finemap_susie_chr21_ENSG00000142166.13" ../polyfun/finemap_susie_rss.sh 21 34567305 ENSG00000142166.13
# sbatch --job-name="finemap_susie_chr21_ENSG00000185917.14" ../polyfun/finemap_susie_rss.sh 21 34567305 ENSG00000185917.14





tail -n +2 /gpfs/commons/home/tlin/script/coloc/snp_of_interest.tsv | while IFS=$'\t' read -r SNP CHR BP source
do
  # Run sbatch with the extracted arguments
  sbatch ../polyfun/finemap_susie_rss.sh "$CHR" "$BP" "$source"
done
