import pandas as pd

print("***************************************")
print("start to merge file........")
print("***************************************")

for i in range(2,3):
	
	#import files	

	brain_atac = pd.read_csv("~/polyfun/data/data_backup/brain_atac_polyfun_backup/brain_atac_seq_chr%s.annot.gz"%i, compression= "gzip", sep = "\t")
	brain_atac_ldscores = pd.read_parquet("~/polyfun/data/data_backup/brain_atac_polyfun_backup/brain_atac_seq_chr%s.l2.ldscore.parquet"%i) 
	bl_ldscores = pd.read_parquet("~/polyfun/data/data_backup/baseline_backup/baselineLF2.2.UKB.%s.l2.ldscore.parquet"%i) 
	bl = pd.read_parquet("~/polyfun/data/data_backup/baseline_backup/baselineLF2.2.UKB.%s.annot.parquet"%i)

	# keep weight file the same
 
	# merge microglia data from brain_atac to BaseLine file
	bl['microglia_atac'] = brain_atac['microglia_atac']  
	bl_ldscores['microglia_atac'] = brain_atac_ldscores["microglia_atac"]
	
	# check if this is merged correctly
	try:
		bl['microglia_atac']
	except KeyError:
		print("not merged correctly, please try to merge again")

	# save the merged file to parquet file

	#bl.to_parquet("~/polyfun/data/BL_brain_atac/BL_microglia.chr%s.annot.parquet"%i) 
	#bl_ldscores.to_parquet("~/polyfun/data/BL_brain_atac/BL_microglia.chr%s.l2.ldscore.parquet"%i)  

        #check
        print("check larger for float32")	
	print(np.where(brain_atac['microglia_atac'].values >= np.finfo(np.float32).max))
	print(np.where(brain_atac_ldscores['microglia_atac'].values >= np.finfo(np.float32).max))
        

        print("\n")
        print("check finite....")
        print(np.isfinite(brain_atac['microglia_atac'].any())
        print(np.isfinite(brain_atac_ldscores['microglia_atac'].any())
                	
        	

print("hoory! finish:)")



