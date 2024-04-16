path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/pip_filt'

sumstat='bellenguez'


for thres in 0.8 0.7 0.6 0.5 0.4 0.3 0.2
do
    echo try thres = $thres
    for anno in susie baseline omics_dl omics
    do  
        touch ${output_path}/${sumstat}_${anno}_pip_${thres}.snp
        for chr in {1..22}
        do  
            echo start chr$chr ...
            if [ "$anno" != "susie" ]; then 
                zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk  'NR==FNR {       
                    elements[$0]  # Store elements in an array
                    next
                }($14 in elements) {print $2 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -
            else
                echo run susie
                zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk  'NR==FNR {       
                    elements[$0]  # Store elements in an array
                    next
                }($13 in elements) {print $1 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -

                #cat ${output_path}/${sumstat}_${anno}_pip_${thres}.snp

            fi
        done

        snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_${thres}.snp )
        echo $snp "SNP(s) were included in " $anno
    done    
done



    #echo wrote ${output_path}/${sumstat}_${anno}_pip_${thres}.snp
    

    # zcat $path/${sumstat}_${anno}_chr1.txt.gz| head -1 > ${output_path}/${sumstat}_${anno}_pip_${thres}_analysis.tsv
    # #touch ${output_path}/${sumstat}_${anno}_pip_${thres}.snp_analysis
    # for chr in {1..22}
    # do  
    #     echo start chr$chr ...
    #     zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk 'NR==FNR {         # For the first files
    #         elements[$0]  # Store elements in an array
    #         next
    #     }($14 in elements) {print $0 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}_analysis.tsv"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -

    # done

    # snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_${thres}_analysis.tsv )
    # echo $snp "SNP(s) were included in " $anno
#done

#baseline omics omics_dl
# 0.7 0.6 0.5 0.4 0.3 0.2 omics omics_dl