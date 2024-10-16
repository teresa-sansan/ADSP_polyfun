#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
#output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/pip_filt'
#output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap/remove_index0'
path=$1 ## aggregated finemap result
output_path=$2 ## output
sumstat='bellenguez'
anno=$3

## get snp count for each thres
## extract snp in remove_index 0
#baseline omics_dl omics

for thres in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 
do  
    echo run pip thres = $thres
    touch ${output_path}/${sumstat}_pip_count.txt

    if [ "$anno" != "susie" ]; then 
        cat $output_path/remove_index0/bellenguez_${anno}.txt |awk  'NR==FNR {       
            elements[$0]  # Store elements in an array
            next
        }($17 in elements) {print $2 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -
    
    else
        cat $output_path/remove_index0/bellenguez_${anno}.txt |awk  'NR==FNR {       
                elements[$0]  # Store elements in an array
                next
            }($15 in elements) {print $3 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -
    fi
    snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_${thres}.snp )
    echo $snp "SNP(s) were included in" $anno, pip threshold =  $thres
    echo $anno $thres $snp >> ${output_path}/${sumstat}_pip_count.txt
    #snp_count=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_count.txt )

#echo $snp_count $thres >> ${output_path}/${sumstat}_${anno}_pip_count.txt

done    



# for anno in susie baseline omics_dl  omics
# do
#     for thres in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 
#     do  
#         echo try thres = $thres
#         touch ${output_path}/${sumstat}_${anno}_pip_count.txt

#         if [ "$anno" != "susie" ]; then 
#             cat ${output_path}/bellenguez_${anno}.txt| awk -v pip=$thres '($11 > pip) {print $2}' > ${output_path}/${sumstat}_${anno}_pip_count.txt1
            
#         else
#             echo run susie 
#             cat ${output_path}/bellenguez_${anno}.txt| awk -v pip=$thres '($10 > pip) {print $1}' > ${output_path}/${sumstat}_${anno}_pip_count.txt1
#         fi
#         snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_count.txt1 )

#     echo $snp $thres >> ${output_path}/${sumstat}_${anno}_pip_count.txt
#     done    


#     rm ${output_path}/${sumstat}_${anno}_pip_count.txt1
#     cat ${output_path}/${sumstat}_${anno}_pip_count.txt
# done



# ## extract snp in pip_thres

# for thres in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 
# do
#     echo try thres = $thres
#     for anno in susie baseline omics_dl omics
#     do  
#         touch ${output_path}/${sumstat}_${anno}_pip_${thres}.snp
#         for chr in {1..22}
#         do  
#             echo start chr$chr ...
#             if [ "$anno" != "susie" ]; then 
#                 zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk  'NR==FNR {       
#                     elements[$0]  # Store elements in an array
#                     next
#                 }($17 in elements) {print $2 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -
#             else
#                 echo run susie
#                 zcat $path/${sumstat}_${anno}_chr${chr}.txt.gz |awk  'NR==FNR {       
#                     elements[$0]  # Store elements in an array
#                     next
#                 }($15 in elements) {print $1 >> "'"${output_path}/${sumstat}_${anno}_pip_${thres}.snp"'" } ' ${output_path}/${sumstat}_${anno}_pip_thres${thres}.txt -

#                 #cat ${output_path}/${sumstat}_${anno}_pip_${thres}.snp

#             fi
#         done

#         snp=$(wc -l < ${output_path}/${sumstat}_${anno}_pip_${thres}.snp )
#         echo $snp "SNP(s) were included in " $anno
#     done    
# done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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