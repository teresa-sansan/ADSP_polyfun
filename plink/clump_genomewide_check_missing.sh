cd bellenguez
touch missing.txt
for dir in ADSP ADSP_qc_all ADSP_qc_variant ADSP_UKBB ADSP_UKBB_qc
#for dir in ADSP
do
 echo $dir| tee -a  missing.txt
 for i in {1..22}
 do
  if [ ! -f "$dir/${dir}_all_qc_${i}.clumped" ]; then
     if [ ! -f "$dir/${dir}_qc_${i}.clumped" ]; then
	if [ ! -f "$dir/${dir}_only_qc_${i}.clumped" ]; then
          echo $dir/${dir}_all_qc_${i} | tee -a missing.txt ; 
        fi
     fi
  fi
 done
 echo >> missing.txt
done

