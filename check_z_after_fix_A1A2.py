import pandas as pd

kunkle_org = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_etal_Stage1.parquet')
kunkle_fix = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/kunkle_2019/AD_Kunkle_polyfun_fixed.parquet')

def check_z(BP):
	chr=kunkle_org[kunkle_org['BP'] == BP]['CHR'].values[0]
	UKBB = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.%s.annot.parquet"%chr)
	print(shape(UKBB))
	print("kunkle_original")
	print(kunkle_org[kunkle_org['BP'] == BP])

	print("\n")
	print("kunkle_fix")
	print(kunkle_fix[kunkle_fix["BP"] == BP])
	print("\n")
	print("UKBB")
	print(UKBB8[UKBB8["BP"] == BP].iloc[:,0:5])



