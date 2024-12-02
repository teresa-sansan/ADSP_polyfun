module purge
module load rstudio-server/2022.12.0-353
module load gcc/11.2.0
module load clonemap
echo "connecting R-4.2.1 port 8777" 
rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.2.1 --www-port 8788 --server-data-dir /nfs/scratch/rstudio/tlin --database-config-file /nfs/sw/rstudio-server/rstudio-server-2022.12.0-353/etc/rstudio/database.conf  --server-user tlin
#rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.2.1 --www-port 1333 --server-data-dir /nfs/scratch/rstudio/tlin --database-config-file /gpfs/commons/home/thasan/easybuild/software/RStudio-Server/2022.07.2+576-foss-2022a-Java-11-R-4.2.1/database.conf  --server-user tlin

                         