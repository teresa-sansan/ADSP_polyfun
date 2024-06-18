import pandas as pd
from functools import reduce  
from tqdm import tqdm


bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.'
deepsea='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/deepsea/deepsea_high_h2_chr'
roadmap='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/roadmap/roadmap_high_h2_chr'
glass_lab_enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer/glass_lab_enformer_high_h2_chr'
glass_lab_enformer_abs='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab_enformer_abs/glass_lab_enformer_abs_high_h2_chr'
enformer='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/enformer/enformer_high_h2_chr'
glass_lab='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/annotations_high_h2/glass_lab/glass_lab_high_h2_chr'


suffix= '.annot.gz'
output_path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC/chr'

def merge_annotation(chrom):
    annotations = [bl, deepsea, roadmap, glass_lab_enformer, glass_lab_enformer_abs, enformer, glass_lab]
    annotation_chr=[]
    for i in range(0, len(annotations)):
        if annotations[i] == '/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/baselineLF2.2.UKB.':
            annotation_name = annotations[i] + str(chrom) + '.annot.parquet'
            anno = pd.read_parquet(annotation_name, engine='pyarrow')
        else:
            annotation_name = annotations[i] + str(chrom) + suffix
            anno = pd.read_csv(annotation_name, compression='gzip', sep = '\t')
        
        
        anno = anno.drop(anno.columns[[0,2,3,4]], axis = 1)
        annotation_chr.append(anno)
            
    anno_merge = reduce(lambda left, right:pd.merge(left,right,on=["SNP"]),annotation_chr)
    anno_merge.insert(1, 'Intercept', 1)
    print(anno_merge.shape)
    anno_merge.to_csv(output_path+str(chrom)+".tsv",index = False,sep='\t')
   
        

#for i in tqdm(range(1,22)):
#    print('running on chr%s'%i)
merge_annotation(9)
    


