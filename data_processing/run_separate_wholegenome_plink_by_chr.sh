#!/bin/bash

for bfile in ADSP/ADSP_all ADSP_qc_all/ADSP_qc_all ADSP_qc_variant/ADSP_qc_variant ADSP_UKBB/ADSP_UKBB_only ADSP_UKBB_qc/ADSP_UKBB_qc_all
do
sbatch --export=bfile=$bfile separate_wholegenome_plink_by_chr.sh
done
