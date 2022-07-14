## Teresa Lin 
## Last edit: 07.11.2022
## Purpose: Replace the previous not converging region with the one that's converge. 

import pandas as pd
import numpy as np


## variables to set before start
fix_convergence_path = '/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/try_rescue_not_converge/aggregate.all.txt.gz'
agg_file_before_convergence_fix = '/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/'
save_file = '/gpfs/commons/home/tlin/output/wightman/fixed_0224/finemap/max_snp_10/agg_fixed_converge.tsv.gz'

## Read file
fix_convergence = pd.read_csv(fix_convergence_path, compression = 'gzip', sep = '\t')
converge_genomewide=pd.DataFrame()  ## create an empty df so it can be merged later.

for i in range(1, 23):
    # print("start chr %s"%(i))
    ori_filename = agg_file_before_convergence_fix + chr{}.aggregrate.all.txt.gz'  ## should be "aggregate" but i mis-spelled. Oops
    ori = pd.read_csv(ori_filename.format(i), sep = '\t')
    converge_chr=ori[~ori.SNP.isin(fix_convergence.SNP)]
    converge_genomewide = pd.concat([converge_genomewide, converge_chr])
    

finemap_genomewide = pd.concat([converge_genomewide, fix_convergence])
finemap_genomewide.CHR = finemap_genomewide.CHR.astype(int)
finemap_genomewide = finemap_genomewide.sort_values("CHR")


duplicated = finemap_genomewide[finemap_genomewide.SNP.duplicated()].SNP
assert len(duplicated) == 0, "There are duplicated SNPs when merging."   ## check theres nothing wrong when running the process. 


finemap_genomewide["POS"] = 'Chr'+ finemap_genomewide.CHR.astype(str) + '_' + finemap_genomewide.start.astype(str) + '_'+ finemap_genomewide.end.astype(str)
finemap_genomewide.to_csv(save_file, sep = '\t', compression = 'gzip',index = False)
print("saved aggregated file in ", save_file)
