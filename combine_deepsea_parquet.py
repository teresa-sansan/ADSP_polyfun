import pandas as pd
import numpy as np
import functools

for k in range(1,23,1):
	print('start at chr%d ...'%k)
	print('start reading files')
	DNase = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr%d.l2.ldscore.parquet'%k)
	H3K27ac = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K27ac/roadmap_H3K27ac_chr%d.l2.ldscore.parquet'%k)
	H3K4me1 = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me1/roadmap_H3K4me1_chr%d.l2.ldscore.parquet'%k)
	H3K4me3 = pd.read_parquet('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me3/roadmap_H3K4me3_chr%d.l2.ldscore.parquet'%k)

	DNase_M = np.loadtxt('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr%d.l2.M'%k, dtype = str)
	H3K27ac_M = np.loadtxt('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K27ac/roadmap_H3K27ac_chr%d.l2.M'%k, dtype = str)
	H3K4me3_M = np.loadtxt('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me3/roadmap_H3K4me3_chr%d.l2.M'%k, dtype = str)
	H3K4me1_M = np.loadtxt('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me1/roadmap_H3K4me1_chr%d.l2.M'%k, dtype = str)


	anno_all = [DNase[['CHR','SNP','BP','A1', 'A2', 'E124-DNase']], 
			H3K27ac[['CHR','SNP','BP','A1','A2','E029-H3K27ac']],
			H3K4me1[['CHR','SNP','BP','A1','A2','E029-H3K4me1']],
			H3K4me3[['CHR','SNP','BP','A1','A2','E029-H3K4me3']]]

	## l2.ldscore.parquet file
	#l2_score = functools.reduce(lambda left,right: pd.merge(left,right, on = ["SNP","CHR","BP","A1","A2"],how = 'outer'),anno_all)
	
	## l2.M
	anno = ['E124-DNase','E029-H3K27ac','E029-H3K4me1','E029-H3K4me3']
	
	ld_M_return = np.array([])
	ld = [DNase, H3K27ac, H3K4me1, H3K4me3]

	ld_M = [DNase_M, H3K27ac_M, H3K4me1_M, H3K4me3_M]
	
	for score, M, anno_index in zip(ld, ld_M, anno):
		for i in range(len(score.columns)):
			if (score.columns[i] == anno_index):
				ld_M_return = np.append(ld_M_return, M[i-5])
				#print(ld_M_return)

	ld_M_return = np.reshape(ld_M_return,(1,6))

	print("start saving file")
	
	np.savetxt('data/roadmap/roadmap_all_anno.%d.l2.M'%k, ld_M_return, delimiter = '\t',fmt='%s')	
	l2_score.to_parquet('~/polyfun/data/roadmap/roadmap_all_anno.%d.l2.ldscore.parquet'%k)
