cd /gpfs/commons/home/tlin/output/sbayesR/fixed_0224
rm bellenguez_agg_snpRes
cat bellenguez_chr1.snpRes| head -1 > bellenguez_agg_snpRes
cat bellenguez_chr*.snpRes |grep -v 'A1'| tr -s ' '| cut -d ' ' -f 2-12  >> bellenguez_agg_snpRes
cat bellenguez_agg_snpRes|sed '1c Id SNP CHR BP A1 A2 A1FREQ BETA SE PIP LastSampleEff' > file_changetitle && mv file_changetitle bellenguez_agg_snpRes

