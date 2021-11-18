zcat finemap_UKBBbaseline.extract_e-01.gz| awk '$9<1e-7 && $10>0.1 {print $0}' > kunkle_rigorous.tsv

