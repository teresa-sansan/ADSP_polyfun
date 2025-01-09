import pandas as pd
#bl_path='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB'
#file = input()
#print('baselineLF2.2.UKB.%d.annot.parquet'%(22))

for i in range(1,22):
    print('reading chr %d...'%(i))
    df = pd.read_parquet('baselineLF2.2.UKB.%d.annot.parquet'%(i))
    
    df.to_csv('/gpfs/commons/home/tlin/data/ukbb_anno/baselineLF2.2.UKB.%d.annot.tsv'%(i), index = False, sep = '\t')
    print(df.shape)
#print('done')


