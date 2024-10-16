import pandas as pd
import numpy as np

## note (08.24.22)
## the P should be MAF instead!
## see ~/script/data_processing/wightman_MAF_beta.ipynb for more detail.
## this is wrong

gwas = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/Wightman_et_al_2021_hg37_ldsc.tsv.gz", sep = '\t', compression = 'gzip')
gwas["BETA"] = gwas.Z/np.sqrt(2*gwas.P*(1-gwas.P)*gwas.N + gwas.Z **2 )
gwas.to_csv('processed/Wightman_2021_hg37_withbeta.tsv', sep = '\t', index = False)
