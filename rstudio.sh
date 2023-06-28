module purge
module load rstudio-server/1.4.1717
module load gcc/11.2.0
module load clonemap
echo "connecting R-4.1.3 port 8777"
rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.1.3  --www-port 8777 --server-data-dir /nfs/scratch/rstudio --database-config-file /nfs/sw/rstudio-server/rstudio-server-1.4.1717/etc/rstudio/database.conf
