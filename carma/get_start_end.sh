#!/bin/bash
cd /gpfs/commons/home/tlin/output/CARMA/geno_filt
# Define range for chr
for chr in {1..22}; do
   echo running chr$chr
    for file in ${chr}_*.gz; do
        if [[ -f $file ]]; then
            ld=$(echo "$file" | cut -d'_' -f2 | cut -d'.' -f1)
        
            start=$(zcat "$file" | tail -n+2 | head -1 | cut -f 3)
            end=$(zcat "$file" | tail -n-1 | cut -f 3)
            
            # Print the output in tab-separated format
            echo -e "${chr}\t${ld}\t${start}\t${end}" >> ../chr_start_end.tsv 
        else
            echo "No files found for chr ${chr}."
        fi
    done
done
