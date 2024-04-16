import pandas as pd
import os
from functools import reduce  


path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/36K_QC/annotated_hg37_plink_ibd/plink_prs/clump_wightman/'
# for chrom in range(5,23):
#     print('start chr%s'%chrom)
#     file_name=path +'ADSP_qc_plink_pT_%s.0.5.profile'%chrom
#     print(file_name)
#     # no_dup=file_name+"_no_dup_space.prs"
#     script="cat "+ file_name + "|tr -s ' '| cut -d ' ' -f 3,7 > %schr%s.score"%(path,chrom)
#     os.system(script)


# Initialize an empty list to store DataFrames
dfs = []

for chrom in range(1,23):
        print('start chr%s'%chrom)
        filename = "chr%d.score"%chrom
        filepath = os.path.join(path, filename)
        df = pd.read_csv(filepath, sep = ' ', names = ['IID','SCORE'])
        dfs.append(df)
        

# Merge the DataFrames on the 'IID' column
prs_merge = reduce(lambda left, right:pd.merge(left,right,on=["IID"]),dfs)
prs_merge = prs_merge.iloc[1:]
prs_merge.columns = ['IID'] + [f'SCORE_{chrom}' for chrom in range(1, 23)]


## merged_with whole genome
prs_wholegenome = pd.read_csv(path + 'pT_0.5.prs', sep = ' ', names = ['IID','SCORE_wholegenome'])
prs_merge = pd.merge(prs_wholegenome,prs_merge, on='IID')

pheno=pd.read_csv('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_ADSP_IBD.tsv', sep='\t')
pheno_merge = pd.merge(pheno,prs_merge, left_on ="SampleID", right_on='IID')


pheno_merge.to_csv(path+ 'pT0.5_all_chrom.tsv', index=False, sep = '\t')


print(pheno_merge)


