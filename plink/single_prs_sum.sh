cd /gpfs/commons/home/tlin/output/prs/bellenguez/updateRSID/interested_SNP/minPIP
file_name='updateRSID_interested_SNP_minPIP_chr'

output_name='updateRSID_interested_SNP_minPIP'
echo write interested.prs

awk '{ sum[$2]+=$6 } END {for (user in sum) print user, sum[user] }' $file_name*.profile > $output_name.prs



