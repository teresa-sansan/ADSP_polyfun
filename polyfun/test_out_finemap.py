## IBSS did not converage
## still need to fix this one
python ~/polyfun_omer_repo/finemapper.py --ld /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr5_102000001_105000001 \
  --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.5.snpvar_constrained.gz \
  --n 487511 \
  --chr 5 --start 102000001 --end 105000001 --method susie \
  --max-num-causal 5 \
  --allow-missing --out test_IBSS_converage.gz




## testing assertion error

## did not happen after using the lastest finemapper.py
python ~/polyfun_omer_repo/finemapper.py --ld /gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/ukb_ld/chr5_102000001_105000001 \
  --sumstats /gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.10.snpvar_constrained.gz \
  --n 487511 \
  --chr 10 --start 103000001 --end 106000001 --method susie \
  --max-num-causal 5 \
  --allow-missing --out test_assertion_error.gz 

