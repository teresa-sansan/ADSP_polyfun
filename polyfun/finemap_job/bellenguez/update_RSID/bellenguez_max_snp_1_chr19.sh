#!/bin/bash 
#SBATCH --job-name=finemap 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=tlin@nygenome.org 
#SBATCH --mem=150G 
#SBATCH --time=70:00:00 
#SBATCH --output=/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/%x_%j.log 

cd ~/polyfun_omer_repo 
source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate 
conda activate polyfun
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 1 --end 3000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.1_3000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 1000001 --end 4000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.1000001_4000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 2000001 --end 5000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.2000001_5000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 3000001 --end 6000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.3000001_6000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 4000001 --end 7000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.4000001_7000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 5000001 --end 8000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.5000001_8000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 6000001 --end 9000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.6000001_9000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 7000001 --end 10000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.7000001_10000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 8000001 --end 11000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.8000001_11000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 9000001 --end 12000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.9000001_12000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 10000001 --end 13000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.10000001_13000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 11000001 --end 14000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.11000001_14000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 12000001 --end 15000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.12000001_15000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 13000001 --end 16000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.13000001_16000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 14000001 --end 17000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.14000001_17000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 15000001 --end 18000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.15000001_18000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 16000001 --end 19000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.16000001_19000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 17000001 --end 20000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.17000001_20000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 18000001 --end 21000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.18000001_21000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 19000001 --end 22000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.19000001_22000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 20000001 --end 23000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.20000001_23000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 21000001 --end 24000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.21000001_24000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 22000001 --end 25000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.22000001_25000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 23000001 --end 26000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.23000001_26000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 24000001 --end 27000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.24000001_27000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 25000001 --end 28000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.25000001_28000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 26000001 --end 29000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.26000001_29000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 27000001 --end 30000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.27000001_30000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 28000001 --end 31000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.28000001_31000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 29000001 --end 32000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.29000001_32000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 30000001 --end 33000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.30000001_33000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 31000001 --end 34000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.31000001_34000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 32000001 --end 35000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.32000001_35000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 33000001 --end 36000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.33000001_36000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 34000001 --end 37000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.34000001_37000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 35000001 --end 38000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.35000001_38000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 36000001 --end 39000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.36000001_39000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 37000001 --end 40000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.37000001_40000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 38000001 --end 41000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.38000001_41000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 39000001 --end 42000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.39000001_42000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 40000001 --end 43000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.40000001_43000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 41000001 --end 44000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.41000001_44000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 42000001 --end 45000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.42000001_45000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 43000001 --end 46000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.43000001_46000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 44000001 --end 47000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.44000001_47000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 45000001 --end 48000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.45000001_48000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 46000001 --end 49000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.46000001_49000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 47000001 --end 50000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.47000001_50000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 48000001 --end 51000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.48000001_51000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 49000001 --end 52000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.49000001_52000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 50000001 --end 53000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.50000001_53000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 51000001 --end 54000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.51000001_54000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 52000001 --end 55000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.52000001_55000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 53000001 --end 56000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.53000001_56000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 54000001 --end 57000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.54000001_57000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 55000001 --end 58000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.55000001_58000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 56000001 --end 59000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.56000001_59000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 57000001 --end 60000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.57000001_60000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 58000001 --end 61000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.58000001_61000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1
python3 /gpfs/commons/home/tlin/polyfun_omer_repo/finemapper_ori.py --chr 19 --start 59000001 --end 62000001 --out /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap_check_pipeline/bellenguez.chr19.59000001_62000001.gz --method susie --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.19.snpvar_constrained.gz --n 487511 --memory 1 --max-num-causal 1