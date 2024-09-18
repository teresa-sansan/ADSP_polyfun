### load all required modules
module purge
module load R/4.2.2
module load gcc/11.2.0
#LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so R
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so 

which R

#/nfs/sw/R/R-4.2.2/bin/Rscript carma.r 22 1
#module load rstudio-server/1.4.1717
#rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.2.2  --www-port 8777 --server-data-dir /nfs/scratch/rstudio --database-config-file /nfs/sw/rstudio-server/rstudio-server-1.4.1717/etc/rstudio/database.conf

module load rstudio-server/2022.12.0-353
echo launch R 4.2.2 on port 8899
rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.2.2 --www-port 8899 --server-data-dir /nfs/scratch/rstudio/tlin --database-config-file /nfs/sw/rstudio-server/rstudio-server-2022.12.0-353/etc/rstudio/database.conf  --server-user tlin

