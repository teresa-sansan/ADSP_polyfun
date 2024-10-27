import os
import sys
import pandas as pd


def process_parquet(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"The file {file} does not exist.")
    
    df = pd.read_parquet(file)
    bp_mean = df['BP'].mean().astype(int)
    chr=df.CHR.values[1]  
    return chr,bp_mean 

if __name__ == "__main__":
    gene = sys.argv[1]  
    file='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene/'+gene+'_sumstats.hg38.parquet'
    chr,mean_value = process_parquet(file)
    print(chr,mean_value)
