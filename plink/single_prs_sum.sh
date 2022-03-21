cd /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP

file_name='updateRSID_interested_SNP_chr'
output_name='updateRSID_interested_SNP'
echo write interested.prs

awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $file_name*.profile > $output_name.prs



