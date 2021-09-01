import pandas as pd
import functools

for i in range(1,23,1):
	print('start at chr%d ...'%i)
	print('start reading files')
	DNase = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/DNASE/expecto_DNASE_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K27ac = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/H3K27ac/expecto_H3K27ac_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K27me3 = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/H3K27me3/expecto_H3K27me3_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K4me1 = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/H3K4me1/expecto_H3K4me1_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K4me3 = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/H3K4me3/expecto_H3K4me3_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)
	H3K9ac = pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/expecto_annotations/annotations_combined_ukb/split/H3K9ac/expecto_H3K9ac_%d.annot.gz'%i,
		compression ='gzip', sep = '\t', header = 0)


	anno_all = [DNase[['CHR','SNP','BP','A1', 'A2', '1932']], 
			H3K27ac[['CHR','SNP','BP','A1','A2','1165']],
			H3K27me3[['CHR','SNP','BP','A1','A2','1656']],
			H3K4me1[['CHR','SNP','BP','A1','A2','1168']],
			H3K4me3[['CHR','SNP','BP','A1','A2','1939']],
			H3K9ac[['CHR','SNP','BP','A1','A2','1941']] ]


	data = functools.reduce(lambda left,right: pd.merge(left,right, on = ["SNP","CHR","BP","A1","A2"],how = 'outer'),anno_all)
	print("start saving file")

	data.to_csv('~/polyfun/data/deepsea/deepsea_all_anno.%d.annot.gz'%i, compression = 'gzip', sep = '\t', index = False)
