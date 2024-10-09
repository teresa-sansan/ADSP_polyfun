import pandas as pd

ld = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/All_Regions_ALL_ensemble_1000G_hg38_EUR.bed', sep = '\t')
ld = ld.rename(columns={'chr':'Chrom','start':'StartBP','stop':'EndBP'})
ld.Chrom = ld.Chrom.str.replace('chr', '')
ld['Block'] = ld.index+1
print(ld.loc[:,['Block','Chrom','StartBP','EndBP']])
ld.loc[:,['Block','Chrom','StartBP','EndBP']].to_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/ref.pos', sep = '\t', index = False)