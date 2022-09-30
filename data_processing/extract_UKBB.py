import pandas as pd
f = open('/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/UKB_rsid.txt',"w")
for chr in range(22):
	chr = chr+1
	print("writing chr " + str(chr))
	file = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB."+str(chr) + ".annot.parquet")
	if( chr==1 ):
		file.iloc[:,:5].to_csv('/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/UKB_rsid.txt' , mode = "a",index = False, sep = '\t',header = True)
	else:
		file.iloc[:,:5].to_csv('/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/UKB_rsid.txt' , mode = "a",index = False, sep = '\t',header = False)

f.close()

## to clear out those that dont follow the rsid rules
!cat UKB_rsid.txt | awk '$2 ~/rs/' > UKB_rsid_clean.txt