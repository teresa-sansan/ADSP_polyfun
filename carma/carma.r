### load all required modules
module purge
module load R/4.1.1
module load gcc/11.2.0
LD_PRELOAD=/gpfs/commons/home/tlin/miniconda3/envs/carma/lib/libmkl_rt.so R

which R
