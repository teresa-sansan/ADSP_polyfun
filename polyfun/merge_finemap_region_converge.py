## Teresa Lin 
## Last edit: 08.15.2022
## Purpose: Replace the regions that encountered IBSS not converging warning with the one that's rescued

import pandas as pd
import numpy as np

## set PATH before start
if False:
    fix_convergence_path = '/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl_brain_atac/max_snp_10/'
    agg_file_before_convergence_fix = '/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl_brain_atac/max_snp_10/'
    save_file = '/gpfs/commons/home/tlin/output/bellenguez/old/bellenguez_fixed_0224_annotations/bl_brain_atac/max_snp_10/agg_fixed_converge.tsv.gz'

if False:
    path = '/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/bl/finemap/max_snp_10/'

if True:
    path='/gpfs/commons/home/tlin/output/jansen/finemap/max_snp_10/'

## Read file
fix_convergence = pd.read_csv((path+'try_rescue_not_converge/aggregate_rescue.all.txt.gz'), compression = 'gzip', sep = '\t')
converge_genomewide=pd.DataFrame()  ## create an empty df so it can be merged later.

## loop over chr1-22
for i in range(1, 23):  
    print("start chr %s"%(i))
    ori_filename = path + 'chr%s.aggregate.all.txt.gz'%(i) 
    ori = pd.read_csv(ori_filename.format(i), sep = '\t')
    converge_chr=ori[~ori.SNP.isin(fix_convergence.SNP)]
    converge_genomewide = pd.concat([converge_genomewide, converge_chr])
    

print('finished gettin the regions were found in both ori and fixed_converging file, now keeping only the one in fixed converge file....')
finemap_genomewide = pd.concat([converge_genomewide, fix_convergence])
finemap_genomewide.CHR = finemap_genomewide.CHR.astype(int)
finemap_genomewide = finemap_genomewide.sort_values("CHR")

print('double check there arent duplicated SNP in the agg file')
duplicated = finemap_genomewide[finemap_genomewide.SNP.duplicated()].SNP
assert len(duplicated) == 0, "There are duplicated SNPs when merging."   ## check theres nothing wrong when running the process. 

finemap_genomewide["POS"] = 'Chr'+ finemap_genomewide.CHR.astype(str) + '_' + finemap_genomewide.start.astype(str) + '_'+ finemap_genomewide.end.astype(str)
save_file = path+'agg_fixed_converge.tsv.gz'
finemap_genomewide.to_csv(save_file, sep = '\t', compression = 'gzip',index = False)
print("saved aggregated file in ", save_file)
