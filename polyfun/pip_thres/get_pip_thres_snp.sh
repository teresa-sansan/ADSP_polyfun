## will create a count file with $anno $pip_thres #snp_that_passed_thres #snp_in_the_credibleset_that_passed_thres

path=$1
anno=$2
#path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa_update'

touch $path/pip_thres/credible_set_snp.count
#baseline omics omics_dl 

echo start $anno
if [ $anno == "susie" ]; then
        pip=9
        credible_set=12
        snp=3
    else
        pip=11
        credible_set=14
        snp=2
    fi
    for thres in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
    do
        ## extract credible set
        zcat $path/bellenguez_remove_index0_${anno}.txt.gz | awk -v pip=$pip -v thres=$thres -v cb=$credible_set -v snp=$snp '$pip > thres {print $snp,$cb}' > $path/pip_thres/${anno}_pip${thres}_credible_set_passed.txt
        count_one_snp=$(cat $path/pip_thres/${anno}_pip${thres}_credible_set_passed.txt| tail -n +2| wc -l| awk '{print $1}') ## count only the snps that passed threshold
        cat $path/pip_thres/${anno}_pip${thres}_credible_set_passed.txt| tail -n+2 |cut -d ' ' -f 2 |sort|uniq > $path/pip_thres/${anno}_pip${thres}_credible_set.txt

        ## extract the snp
        file1="${path}/bellenguez_remove_index0_${anno}.txt.gz"
        file2="${path}/pip_thres/${anno}_pip${thres}_credible_set.txt"

        ## create count
        zcat "$file1" | awk -v cb=$credible_set 'NR==FNR{a[$1];next} ($cb in a)' "$file2" - | cut -f $snp > $path/pip_thres/${anno}_pip${thres}.snp  
        count_all_snp=$(wc -l $path/pip_thres/${anno}_pip${thres}.snp| awk '{print $1}') ## count all snps in the credible set
        echo $anno $thres  $count_one_snp $count_all_snp| tee -a  $path/pip_thres/credible_set_snp.count
    done 
