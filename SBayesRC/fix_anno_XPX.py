import pandas as pd


path='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/4SBayesRC'

non_binary = ['SNP', 'Intercept', 'GERP.NS_lowfreq', 'GERP.NS_common',
       'MAF_Adj_Predicted_Allele_Age_common', 'MAF_Adj_LLD_AFR_lowfreq',
       'MAF_Adj_LLD_AFR_common', 'Recomb_Rate_10kb_lowfreq',
       'Recomb_Rate_10kb_common', 'Nucleotide_Diversity_10kb_lowfreq',
       'Nucleotide_Diversity_10kb_common', 'Backgrd_Selection_Stat_lowfreq',
       'Backgrd_Selection_Stat_common', 'CpG_Content_50kb_lowfreq',
       'CpG_Content_50kb_common', 'MAF_Adj_ASMC_lowfreq',
       'MAF_Adj_ASMC_common', 'GTEx_eQTL_MaxCPP_common',
       'BLUEPRINT_H3K27acQTL_MaxCPP_common',
       'BLUEPRINT_H3K4me1QTL_MaxCPP_common',
       'BLUEPRINT_DNA_methylation_MaxCPP_common',
       'Human_Enhancer_Villar_Species_Enhancer_Count_lowfreq',
       'Human_Enhancer_Villar_Species_Enhancer_Count_common']

for chrom in range(1,20):
    print("start chr%d"%chrom)
    anno = pd.read_csv(path + '/chr'+str(chrom)+'.tsv', sep = '\t')
    col_to_remove = ['Coding_UCSC_lowfreq','Coding_UCSC_common']
    #binary_cols = [col for col in chr21_all.columns if chr21_all[col].nunique() == 2 ]
    anno = anno.drop(columns=col_to_remove, axis=1)
    anno.to_csv(path+'chr%d_fixed.tsv'%chrom, sep = '\t', index = False)


