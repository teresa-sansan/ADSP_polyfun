import pandas as pd
import numpy as np



for i in range(1,22):
		
	#import files	
	brain_atac = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr%s.annot.gz"%i, compression= "gzip", sep = "\t")
	brain_atac_ldscores = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr%s.l2.ldscore.parquet"%i) 
	bl_ldscores = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.%s.l2.ldscore.parquet"%i) 
	bl = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.%s.annot.parquet"%i)

	hist_list = ['DNASE', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K9me3', 'HISTONE_OTHER', 'TF']
	for k in hist_list:
		deepsea_ldscores = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/%s/expecto_%s_%s.l2.ldscore.parquet"%(k,k,i))	
		deepsea  = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/%s/expecto_%s_%s.annot.gz"%(k,k,i), compression = "gzip", sep = "\t")	
	
 
		#deepsea  = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/HISTONE_OTHER/expecto_HISTONE_OTHER_%s.annot.gz"%i, compression = "gzip", sep = "\t")
		#deepsea_ldscores = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/HISTONE_OTHER/expecto_HISTONE_OTHER_%s.l2.ldscore.parquet"%i)	
	

		#roadmap = pd.read_csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr%s.annot.gz"%i, sep= '\t')
		#roadmap_ldscores = pd.read_parquet("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr%s.l2.ldscore.parquet"%i) 	
	
		if(brain_atac.shape[0] != bl.shape[0] or deepsea.shape[0] != bl.shape[0]):
			print("annotation inconsistancy, with")
			print("bl.anno = %d, brain_atac.anno = %d, deepsea.anno = %d" %(bl.shape[0], brain_atac.shape[0], deepsea.shape[0]))
	
		if(brain_atac_ldscores.shape[0] != bl_ldscores.shape[0] or deepsea_ldscores.shape[0] != bl_ldscores.shape[0]):
			print("ldscore inconsistancy, with")
			print('chr %s,'%i)
			print('annotation: %s'%k)
			print("bl.ld = %d, brain_atac.ld = %d, deepsea.ld = %d" %(bl_ldscores.shape[0], brain_atac_ldscores.shape[0], deepsea_ldscores.shape[0]))
			print('\n')
	#if(bl.shape[0] != bl_ldscores.shape[0]):
	#	print("bl.anno = %d, bl.ldscores = %d" %(bl.shape[0], bl_ldscores.shape[0]))
print("hoory! finish:)")


# test the NaN warning in running polyfun 1-2
# ValueError: input contains NaN, infinity or a valuetoo large for dtype('float32')

