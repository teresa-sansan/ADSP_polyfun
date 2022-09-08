import pandas as pd
chr22_anno = pd.read_parquet('baselineLF2.2.UKB.22.annot.parquet')
chr22_ldscore = pd.read_parquet('baselineLF2.2.UKB.22.l2.ldscore.parquet')
MAF = []

for i in chr22_ldscore.columns:
	if ("MAF" in i):
		MAF.append(i)


print(chr22_ldscore.loc[:5,MAF])	
