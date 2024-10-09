#sumstat_file='Bellenguez_et_al_2021_hg37_new_sep20_qc.tsv'

sumstat_file=''
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed
chr=$1
echo start chr$chr
echo "SNP A1 A2 freq b se p N" |tr ' ' '\t'>  sbayesrc/bellenguez_hg19_chr$chr.cojo
cat $sumstat_file| awk -v chr="$chr" '$2 == chr {print $1, $4, $5, $7,$11,$12,$6,$19}' |tr ' ' '\t' >> sbayesrc/bellenguez_hg19_chr$chr.cojo