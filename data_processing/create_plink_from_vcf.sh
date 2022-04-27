cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/17K_final/annotated_filtered_hg37

for i in {1..22}
do
echo running chr $i
sbatch --mem=70G -c 5 --wrap="~/plink --vcf ADSP_annotated_chr$i.vcf.gz --recode --double-id -allow-extra-chr --out plink/ADSP_annotated_chr$i"
done


## added two flags (--double-id & --allow_extra-chr) because
## there were and 
## Error: Invalid chromosome code 'chr17_gl000205_random' on line 1790 of .vcf file. (e.g. in chr21)
## increase the memrory size. ( bus error)
