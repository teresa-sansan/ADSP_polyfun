import pandas as pd
import functools

for i in range(1,23,1):
	print('start at chr%d ...'%i)
	print('start reading files')
	DNase = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K27ac = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K27ac/roadmap_H3K27ac_chr%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K4me1 = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me1/roadmap_H3K4me1_chr%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K4me3 = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me3/roadmap_H3K4me3_chr%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)


	anno_all = [DNase[['CHR','SNP','BP','A1', 'A2', 'E124-DNase']], 
			H3K27ac[['CHR','SNP','BP','A1','A2','E029-H3K27ac']],
			H3K4me1[['CHR','SNP','BP','A1','A2','E029-H3K4me1']],
			H3K4me3[['CHR','SNP','BP','A1','A2','E029-H3K4me3']]]


	data = functools.reduce(lambda left,right: pd.merge(left,right, on = ["SNP","CHR","BP","A1","A2"],how = 'outer'),anno_all)
	print("start saving file")

	data.to_csv('~/polyfun/data/roadmap/roadmap_all_anno.%d.annot.gz'%i, compression = 'gzip', sep = '\t', index = False)
