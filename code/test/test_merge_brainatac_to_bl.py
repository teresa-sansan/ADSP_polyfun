import pandas as pd
import numpy as np

print("***************************************")
print("start debugging........")

floatcheck = [] 
nan = []
infinite = []
for i in range(1,23):
		
	#import files	
	print("start checking chr%s"%i)
	brain_atac = pd.read_csv("~/polyfun/data/data_backup/brain_atac_polyfun_backup/brain_atac_seq_chr%s.annot.gz"%i, compression= "gzip", sep = "\t")
	brain_atac_ldscores = pd.read_parquet("~/polyfun/data/data_backup/brain_atac_polyfun_backup/brain_atac_seq_chr%s.l2.ldscore.parquet"%i) 
	bl_ldscores = pd.read_parquet("~/polyfun/data/data_backup/baseline_backup/baselineLF2.2.UKB.%s.l2.ldscore.parquet"%i) 
	bl = pd.read_parquet("~/polyfun/data/data_backup/baseline_backup/baselineLF2.2.UKB.%s.annot.parquet"%i)

	# keep weight file the same
 	
	if(brain_atac.shape[0] != bl.shape[0]):
		print("orignal fail")
	if(brain_atac_ldscores.shape[0] != bl_ldscores.shape[0]):
		print("ldscore fail")
		print(bl_ldscores.shape[0],brain_atac_ldscores.shape[0])
	if(brain_atac_ldscores.shape[0] != brain_atac.shape[0]):
		print("dif. num of rows in anno and ld files from brain_atac_seq, ")
		print(brain_atac_ldscores.shape[0], brain_atac.shape[0])
	# merge microglia data from brain_atac to BaseLine file
	#bl['microglia_atac'] = brain_atac['microglia_atac']  
	#bl_ldscores['microglia_atac'] = brain_atac_ldscores["microglia_atac"]
	
	# check if this is merged correctly
	#try:
	#	bl['microglia_atac']
	#except KeyError:
	#	print("not merged correctly, please try to merge again")

	# save the merged file to parquet file

	#bl.to_parquet("~/polyfun/data/BL_brain_atac/BL_microglia.chr%s.annot.parquet"%i) 
	#bl_ldscores.to_parquet("~/polyfun/data/BL_brain_atac/BL_microglia.chr%s.l2.ldscore.parquet"%i)  

        ##checking for bugs
	#print("check larger for float32.....")	
	#print(np.where(bl['microglia_atac'].values >= np.finfo(np.float32).max))
	#print(np.where(bl_ldscores['microglia_atac'].values >= np.finfo(np.float32).max))	
		
	#print("check finite....")
	#print(np.isfinite(brain_atac['microglia_atac'].any()))
	#print(np.isfinite(brain_atac_ldscores['microglia_atac'].any()))
	#if(np.isinf(bl['microglia_atac'].any()) == True or np.isinf(bl_ldscores['microglia_atac'].any()) ==True):
		infinite.append("chr%s"%i)
		
	#print("check nan ...")      
	#if (np.any(np.isnan(bl_ldscores['microglia_atac'])) == True or np.any(np.isnan(bl['microglia_atac'])) == True):
		nan.append("chr%s"%i)

#print('too large for float32:')
#print(floatcheck)
		
#print("infinite:")
#print(infinite)

#print("nan:")
#print(nan)
print("hoory! finish:)")



