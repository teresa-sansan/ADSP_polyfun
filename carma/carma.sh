### load all required modules
module purge
module load R/4.1.1
module load gcc/11.2.0
#LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so R
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so


module load rstudio-server/1.4.1717
rserver --secure-cookie-key-file=$HOME/rstudio-server/secure-cookie-key --rsession-which-r /nfs/sw/R/R-4.1.1  --www-port 8777 --server-data-dir /nfs/scratch/rstudio --database-config-file /nfs/sw/rstudio-server/rstudio-server-1.4.1717/etc/rstudio/database.conf
