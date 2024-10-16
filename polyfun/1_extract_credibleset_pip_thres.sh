#This script is to extract the credible set block that has at least one SNP with PIP > threshold

#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap_v2/combined'
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap_v3/SNPVAR_SNPS/'
output='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/PRS/36k_hg38/snp_pip_thres'
sumstat='bellenguez'

for anno in baseline omics_dl omics
do
    echo run $anno
    echo "removing credible set index = 0 ..."
    zcat $path/bellenguez_${anno}_chr22.txt.gz | head -1 > $output/remove_index0/bellenguez_${anno}.txt
    for chr in {1..22}
    do  
        zcat $path/bellenguez_${anno}_chr${chr}.txt.gz | tail -n+2| grep -v '1:0' >> $output/remove_index0/bellenguez_${anno}.txt
    done
    for thres in 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0 
    do     
        touch $output/${sumstat}_${anno}_pip_thres${thres}.txt  
        if [ "$anno" == "susie" ]; then 
            cat $output/remove_index0/bellenguez_${anno}.txt|tail -n+2| awk -v thres="$thres" '$12 > thres {print $15}'| sort | uniq >>  $output/${sumstat}_${anno}_pip_thres${thres}.txt      

        else 
            cat $output/remove_index0/bellenguez_${anno}.txt|tail -n+2| awk -v thres="$thres" '$14 > thres {print $17}'| sort | uniq >>  $output/${sumstat}_${anno}_pip_thres${thres}.txt          
        fi
    done
    bash 2_extract_snp_in_credibleset.sh $path $output $anno
    echo
done

echo finished!



