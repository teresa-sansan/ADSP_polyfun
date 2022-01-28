python re-clone/polypred.py \
  --combine-betas \
  --betas polypred_example/bolt.betas.gz,modified/output/polyfun_output.agg.txt.gz \
  --pheno polypred_example/1000G.pheno.txt \
  --output-prefix modified/tutorial__ \
  --plink-exe ~/plink polypred_example/1000G.subset.*.bed
