path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/pip_filt'

sumstat='bellenguez'
thres=0.25

echo try thres = $thres
for anno in baseline
do  
    touch ${output_path}/${sumstat}_${anno}_pip_${thres}.snptest
    for chr in {1..22}
    do  
        echo start chr$chr ...
        zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk 'NR==FNR {         # For the first file
            elements[$0]  # Store elements in an array
            next
        }($14 in elements) {print $2 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snptest"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -

    done

    #echo wrote ${output_path}/${sumstat}_${anno}_pip_${thres}.snp
    snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_${thres}.snptest )
    echo $snp "SNP(s) were included in " $anno
done

#baseline omics omics_dl
### susie's SNP and CHR are opposite ($1,$2)
### susie need to be $13 instead of $14