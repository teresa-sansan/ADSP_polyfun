#for i in {25391373..25399827}
#for i in {25425368..25479354}
title='Finemap_kunkle'
#for i in {25348771..25566479}
for i in {25479329..25479333}
do

if [ -e ${title}${i}.log ]
then
error=$(cat ${title}${i}.log| grep Assertion| head | wc -l)
 if [ $error != 0 ]
 then
 echo $i, $error 
 else
 echo move $i
 mv ${title}${i}.log log_no_Assertion_error
 fi
fi
done
