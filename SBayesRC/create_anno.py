import pandas as pd
from functools import reduce  
from tqdm import tqdm


#bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.'
bl='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations_hg38/merged_annotations_ADSP_v2/baseline_filtered/baseline_chr'
# deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
# roadmap='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/roadmap/roadmap_high_h2_chr'
# glass_lab_enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer/glass_lab_enformer_high_h2_chr'
# glass_lab_enformer_abs='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer_abs/glass_lab_enformer_abs_high_h2_chr'
# enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/enformer/enformer_high_h2_chr'
# glass_lab='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab/glass_lab_high_h2_chr'

sumstat = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg38_from_parquet.cojo'
suffix= '.annot.gz'
output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/hg38/bl_chr'
#sumstat = '/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/'


# def merge_annotation(chrom):
#     annotations = [bl, deepsea, roadmap, glass_lab_enformer, glass_lab_enformer_abs, enformer, glass_lab]
#     annotation_chr=[]
#     for i in range(0, len(annotations)):
#         if annotations[i] == '/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.':
#             annotation_name = annotations[i] + str(chrom) + '.annot.parquet'
#             anno = pd.read_parquet(annotation_name, engine='pyarrow')
#         else:
#             annotation_name = annotations[i] + str(chrom) + suffix
#             anno = pd.read_csv(annotation_name, compression='gzip', sep = '\t')
        
#         anno['ID'] = anno.apply(lambda row: "%s_%s_%s_%s" % (row['CHR'], row['SNP'], row['A1'], row['A2']), axis=1)
#         anno = anno.drop(anno.columns[[0,1,2,3,4]], axis = 1)
#         #anno = anno.drop(anno.columns[[0,2,3]], axis = 1)
#         annotation_chr.append(anno)
            
#     anno_merge = reduce(lambda left, right:pd.merge(left,right,on=["ID"]),annotation_chr)
#     anno_merge.insert(0, 'Intercept', 1)
#     print(anno_merge.shape)
#     return(anno_merge)
#    # anno_merge.to_csv(output_path+str(chrom)+".tsv",index = False,sep='\t')
   
        

# for i in tqdm(range(1,22)):
#     #print('running on chr%s'%i)
#     cojo = pd.read_csv(sumstat+str(i)+'.cojo', sep = '\t')
#     cojo['ID'] =  cojo.apply(lambda row: "%s_%s_%s_%s" % (i, row['SNP'], row['A1'], row['A2']), axis=1)
#     anno = merge_annotation(i)
#     merged_anno =pd.merge(cojo, anno, on='ID', how='inner')
#     merged_anno = merged_anno.drop(merged_anno.columns[[1,2,3,4,5,6,7,8]], axis = 1)
#     merged_anno.to_csv(output_path+str(i)+".tsv",index = False,sep='\t')


for chrom in tqdm(range(1,23)):
## create bl anno
    annotation_name = bl + str(chrom) + '.annot.gz'
    anno = pd.read_csv(annotation_name, compression='gzip', sep = '\t')
    anno = anno.drop(anno.columns[[0,2,3,4]], axis = 1)
    anno = anno.rename(columns = {"A2":"A1", "A1":"A2"})
    anno.insert(1, 'Intercept', 1)
    anno.to_csv(output_path+str(chrom)+".tsv",index = False,sep='\t')