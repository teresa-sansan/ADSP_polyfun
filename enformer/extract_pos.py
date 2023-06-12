import pandas as pd

gtf=pd.read_csv('/gpfs/commons/home/tlin/data/enformer/gencode.v27.gene.gtf', header=None, sep = '\t')
gtf['pos']=gtf.apply(lambda row: row[3] if row[6] == '+' else row[4], axis=1)
gtf["gene_name"] = gtf[8].str.split(';', expand=True)[2].str.replace('gene_name', '').str.strip().str.replace('"', '')
gtf["gene_type"] = gtf[8].str.split(';', expand=True)[1].str.replace('gene_type', '').str.strip().str.replace('"', '')
gtf["gene_id"] = gtf[8].str.split(';', expand=True)[0].str.replace('gene_id', '').str.strip().str.replace('"', '')
gtf = gtf.rename(columns={0: 'chr', 6: 'direction' })
gtf.iloc[:,[0,6,9,10,11,12]].to_csv('/gpfs/commons/home/tlin/data/enformer/gencode.v27.gene_processed.tsv', sep='\t', index=False)