cd /gpfs/commons/home/tlin/polyfun_omer_repo

source /gpfs/commons/groups/knowles_lab/software/anaconda3/bin/activate
conda activate polyfun

echo "start testing rho (change from 0.95 to 0.5)"
#FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr6_40000001_43000001"
FILES="/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr6_39000001_42000001"
sumstat="/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/bellenguez_all.6.snpvar_constrained.gz"

chr=6

start=39000001
end=42000001
	
python finemapper_rho.py \
	--ld $FILES \
	--sumstats $sumstat \
	--n 487511 \
  	--chr ${chr} --start $start --end $end \
  	--method susie \
	--max-num-causal 10 \
  	--allow-missing \
	--out /gpfs/commons/home/tlin/data/test_rho/finemap_bellenguez.${chr}.$start.$end.gz
