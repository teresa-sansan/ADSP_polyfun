curl -k -H "Content-Type: application/json" -X POST -d '{"snps": "/gpfs/commons/home/tlin/data/clean_bellenguez_ld_snp_only.tsv", "pop": "EUR","r2_d": "d", "genome_build": "grch37"}' https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token=d2e42d726f01


## test if only use first 100 SNPs
curl -k -H "Content-Type: application/json" -X POST -d '{"snps": "/gpfs/commons/home/tlin/data/check_100.tsv", "pop": "EUR","r2_d": "d", "genome_build": "grch37"}' https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token=d2e42d726f01

##output
"error": "Internal server error. Please contact LDlink admin."
