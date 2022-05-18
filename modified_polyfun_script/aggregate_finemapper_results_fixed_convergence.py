## modified the aggregate scrip so that I can aggregate the regions that has 1MB window, 0.5MB overlap instead of 3 MB window and 1 MB overlap
## The idea is the same. I am still taking the most centric one. 

## modified line 45 in the original script.
## changed '%s.chr%s.%s_%s.gz' to the output prefix I set for finemap job

import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import logging
import scipy.stats as stats
from tqdm import tqdm
from polyfun import configure_logger, check_package_versions
from polyfun_utils import set_snpid_index
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from polyfun_utils import DEFAULT_REGIONS_FILE


def main(args):
    
    #read sumstats file
    try:
        df_sumstats = pd.read_parquet(args.sumstats)
    except (ArrowIOError, ArrowInvalid):
        df_sumstats = pd.read_table(args.sumstats, sep='\s+')
        
    #compute p-values if needed
    if args.pvalue_cutoff is not None:
        df_sumstats['P'] = stats.chi2(1).sf(df_sumstats['Z']**2)
        
        
    #read regions file
    #df_regions = pd.read_table(args.regions_file)
    
    table = pd.read_table(args.regions_file, sep ='.', names = ["CHR","POS","else"]) 
    if not os.path.exists(args.regions_file):
       print("regions file not found, please recheck the path.")

    table = table[1:]
    table.index -=1
    new = pd.DataFrame(table.POS.str.split('_').tolist(), columns=["START","END"]).astype(int)
    new.START = new.END - 1000000 ## only take the last 1MB of the 3MB window
    new["CHR"] = table.CHR.str.replace('chr','')
    


    ## caluclate the three start sites of three windows that cover the 1MB window. 

    new['block_1'] = new.START -500000
    new['block_2'] = new.START + 500000   ## a.k.a new.end - 500000 
    new['block_3'] = new.END + 500000 

    new1 = pd.concat([new.block_1, new.block_2, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
    new2 = pd.concat([new.START, new.END, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
    new3 = pd.concat([new.block_2, new.block_3, new.CHR], axis=1,keys=["START","END","CHR"]).sort_values(['CHR',"START"])
    
    
    ## concat all the 1MB window and sort them
    df_regions = pd.concat([new1,new2,new3]) 
    df_regions.astype(int)
    df_regions = df_regions[df_regions.START > 0]
    df_regions = df_regions.sort_values(['CHR',"START"]).reset_index().drop(columns=["index"]) 
#    print(df_regions)

 #   if args.chr is not None:
 #       print(df_regions.head)
 #       df_regions = df_regions.query('CHR==%d'%(args.chr))
 #       if df_regions.shape[0]==0: raise ValueError('no SNPs found in chromosome %d'%(args.chr))

    #df_regions = df_regions.loc[df_regions.apply(lambda r: np.any((df_sumstats['CHR']==r['CHR']) & (df_sumstats['BP'].between(r['START'], r['END']))), axis=1)]

    #aggregate outputs
    df_sumstats_list = []
    logging.info('Aggregating results...')
    for _, r in tqdm(df_regions.iterrows()):
        #chr_num, start, end, url_prefix = r['CHR'], r['START'], r['END'], r['URL_PREFIX']
        chr_num, start, end = r['CHR'], r['START'], r['END']
        r = r.astype(int)
        #apply p-value filter if needed
        if args.pvalue_cutoff is not None:
            df_sumstats_r = df_sumstats.query('CHR==%d & %d <= BP <= %d'%(chr_num, start, end))
            if np.all(df_sumstats_r['P'] > args.pvalue_cutoff): continue        
        #output_file_r = '%s.chr%s.%s_%s.gz'%(args.out_prefix, chr_num, start, end)   ## original script        
        #output_file_r = '%s.%s.%s.%s.gz'%(args.out_prefix, chr_num, start, end)  ##need to change again here
        output_file_r = '%s_chr%s.%s.%s.gz'%(args.out_prefix, chr_num, start, end)   ## original script

	
        if not os.path.exists(output_file_r):
            print('not found %s_chr%s.%s.%s.gz'%(args.out_prefix, chr_num, start, end))
            err_msg = 'output file for chromosome %s bp %s-%s doesn\'t exist'%(chr_num, start, end)
            if args.allow_missing_jobs:
                logging.warning(err_msg)
                logging.warning(output_file_r) ## check path
                continue
            else:
                raise IOError(err_msg + '.\nTo override this error, please provide the flag --allow-missing-jobs')

        df_sumstats_r = pd.read_table(output_file_r)
        #mark distance from center
        middle = (start+end)//2
        df_sumstats_r['DISTANCE_FROM_CENTER'] = np.abs(df_sumstats_r['BP'] - middle)
        df_sumstats_r['start'] = start
        df_sumstats_r['end'] = end
        df_sumstats_list.append(df_sumstats_r)
    if len(df_sumstats_list)==0:
        print(df_sumstats_list)
        raise ValueError('no output files found')
    
    
    #keep only the most central result for each SNP
    df_sumstats = pd.concat(df_sumstats_list, axis=0)
    df_sumstats.sort_values('DISTANCE_FROM_CENTER', inplace=True, ascending=True)
    df_sumstats = set_snpid_index(df_sumstats, allow_duplicates=True)
    df_sumstats = df_sumstats.loc[~df_sumstats.index.duplicated(keep='first')]
    del df_sumstats['DISTANCE_FROM_CENTER']
    df_sumstats.sort_values(['CHR', 'BP'], inplace=True, ascending=True)
    
    #write output file
    if args.adjust_beta_freq:
        df_sumstats['BETA_MEAN'] /= np.sqrt(2*df_sumstats['MAF']*(1-df_sumstats['MAF']))
        df_sumstats['BETA_SD']   /= np.sqrt(2*df_sumstats['MAF']*(1-df_sumstats['MAF']))
    df_sumstats.to_csv(args.out, sep='\t', index=False)
    logging.info('Wrote aggregated results to %s'%(args.out))
        


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    #general parameters
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--out-prefix', required=True, help='prefix of output files')
    parser.add_argument('--out', required=True, help='name of the aggregated output files')
    parser.add_argument('--allow-missing-jobs', default=False, action='store_true', help='whether to allow missing jobs')
    parser.add_argument('--regions-file', default=DEFAULT_REGIONS_FILE, help='name of file of regions and their URLs')
    parser.add_argument('--chr', default=None, type=int, help='Target chromosome (if not provided, all chromosomes will be considered)')
    parser.add_argument('--pvalue-cutoff', type=float, default=None, help='only consider regions that have at least one SNP with a p-value greater than this cutoff')
    parser.add_argument('--adjust-beta-freq', default=False, action='store_true', help='If specified, the posterior estimates of the SNP effect sizes will be on per-allele scale rather than a per-standardized genotype scale')
    
    
    #check package versions
    check_package_versions()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))    

    #configure logger
    configure_logger(args.out_prefix)

    #invoke main function
    main(args)
    
