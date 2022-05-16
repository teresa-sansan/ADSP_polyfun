cd /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/try_rescue_not_converge
touch no_result.txt
touch log_result_for_no_result.txt
while read -r line
do
gzfile=$(echo $line| cut -d '.' -f 1-4)
if [ ! -f $gzfile ];then
#echo $gzfile | tee -a >> no_result.txt
head $line >> log_result_for_no_result.txt 
echo >> log_result_for_no_result.txt

fi
done < logfile.txt
