import pandas as pd
import numpy as np
import sys

chrom = sys.argv[1]
ld_index = sys.argv[2]

print(f'running on chr{chrom}, ld_blk: {ld_index}')

## create index file (the row with more than half nan LD)
path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/'
ld_path = f'/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/chr{chrom}_{ld_index}.ld'

ld = np.genfromtxt(ld_path, delimiter=' ')
nan_count_per_row = np.sum(np.isnan(ld), axis=1)

index = np.where(nan_count_per_row > ld.shape[0]/2)[0]
nan_results = nan_count_per_row[index]

index_df = pd.DataFrame({'Row': index, 'NaN_Count':nan_count_per_row[index]})
output_name = f'/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/count/chr{chrom}_{ld_index}_nan_row.tsv'
index_df.to_csv(output_name, index = False, sep = '\t')
print(f'Created index file with %d SNPs' %index_df.shape[0])
#index_df = pd.read_csv(output_name, sep = '\t')

## merge snp information
snp_file = pd.read_csv(f'{path}/chr{chrom}_{ld_index}.bim', sep='\t', header=None)
snp_maf = pd.read_csv(f'{path}/chr{chrom}_{ld_index}.frq', delim_whitespace=True) ## can't directly use this file because this has all the snp (including the ones we filtered)

snp_file.columns = ['CHR','RSID','CENTIMORGANS','BP','A1','A2']
snp_file.index = np.arange(0, len(snp_file))  
check_snp = snp_file[snp_file.index.isin(index_df.Row)]  
check_snp  = check_snp.join(index_df)

output = pd.merge(check_snp,snp_maf[['SNP','MAF','NCHROBS']], left_on = 'RSID', right_on='SNP') 
output.to_csv(f'{path}count/chr{chrom}_{ld_index}_nan_snp_.tsv', sep = '\t', index = False)