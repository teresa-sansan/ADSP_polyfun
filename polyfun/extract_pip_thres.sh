#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap_v2/combined'
sumstat='bellenguez'

#zcat $path/bellenguez_${anno}_chr${chr}.txt.gz | grep -v '1:0' > $path/remove_index0/bellenguez_${anno}_chr${chr}.txt
zcat $path/bellenguez_${anno}.txt.gz | grep -v '1:0' > $path/remove_index0/bellenguez_${anno}.txt

cd $path
for anno in susie baseline omics_dl omics
do
    for thres in 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
    do
        echo thres = $thres
        touch pip_filt/${sumstat}_${anno}_pip_thres${thres}.txt
        for chr in {1..22}
        do  
            echo starting chr $chr ...
            if [ "$anno" == "susie" ]; then 
                zcat  ${sumstat}_${anno}_chr${chr}.txt.gz |tail -n+2| awk -v thres="$thres" '$10 > thres {print $13}'| sort | uniq >>  pip_filt/${sumstat}_${anno}_pip_thres${thres}.txt
            else 
                zcat  ${sumstat}_${anno}_chr${chr}.txt.gz |tail -n+2| awk -v thres="$thres" '$11 > thres {print $14}'| sort | uniq >>  pip_filt/${sumstat}_${anno}_pip_thres${thres}.txt
            fi

        done
    done
    echo
    echo finished $anno
done


##baseline omics_dl omics susie baseline omics_dl omics

## remember that susie doenst have SNPVAR column 