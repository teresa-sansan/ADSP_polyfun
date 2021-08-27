cat anno_new.txt|while read line
do

if  grep -v -q $line BL_anno.txt  
then 
	echo $line
#elif [ $(grep -c $line BL_anno.txt) -gt 1 ]
elif [ $(grep -o $line BL_anno.txt | wc -l ) -eq 1 ]
then
	echo $line
fi         
done
