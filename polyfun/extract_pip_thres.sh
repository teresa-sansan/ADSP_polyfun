path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap'
sumstat='bellenguez'
cd $path
for anno in baseline omics_dl omics 
do
    touch pip_filt/${sumstat}_${anno}_pip_thres0.1.txt
    for chr in {1..22}
    do  
        echo starting chr $chr ...
        zcat  ${sumstat}_${anno}_chr${chr}.txt.gz |tail -n+2| awk '$11 > 0.25 {print $14}'| sort | uniq >>  pip_filt/${sumstat}_${anno}_pip_thres0.25.txt

    done
    echo
    echo finished $anno
done


##baseline omics_dl omics susie 

## remember that susie doenst have SNPVAR column 