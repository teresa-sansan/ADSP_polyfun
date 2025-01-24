import os
import sys
import pandas as pd


def process_parquet(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"The file {file} does not exist.")
    
    df = pd.read_parquet(file)
    #bp_mean = df['BP'].mean().astype(int)
    #bp_range = df['BP'].max().astype(int) - df['BP'].min().astype(int)
    bp_start = df.BP.values[0]
    chr=df.CHR.values[0]  
    return chr, bp_start
    #return chr,bp_mean 

if __name__ == "__main__":
    gene = sys.argv[1]  
    file='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/eQTL/microglia/gene/'+gene+'_sumstats.hg38.parquet'
    chr,bp_start = process_parquet(file)
    print(chr,bp_start)
