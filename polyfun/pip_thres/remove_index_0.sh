#!/bin/bash
#SBATCH --job-name=bellenguez_index
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa/pip_thres/%x_%j.log

anno=$1
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/finemap_v3_backup_teresa'
# for anno in  susie baseline omics omics_dl
# do  
if [ $anno == "susie" ]; then
    echo running $anno ...
    zcat $path/bellenguez_${anno}.txt.gz  | awk '($12 !~ /1:0/) {print$0 }'| gzip > $path/bellenguez_remove_index0_${anno}.txt.gz
    #zcat $pathbellenguez_susie.txt.gz | grep -v '1:0' | gzip > $path/bellenguez_remove_index0_${anno}.txt.gz
    
    zcat $path/bellenguez_remove_index0_${anno}.txt.gz| wc -l
    
else
    echo running $anno ... 
    zcat $path/bellenguez_${anno}.txt.gz  | awk '($14 !~ /1:0/) {print$0 }'|gzip > $path/bellenguez_remove_index0_${anno}.txt.gz
    
fi
#done