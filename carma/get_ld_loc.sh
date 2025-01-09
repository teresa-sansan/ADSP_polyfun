# chr=$1
# pos=$2

# awk -v chr="$chr" -v pos="$pos" '$1 == chr && pos >= $3 && pos <= $4 { print $1,$2 }' /gpfs/commons/home/tlin/output/CARMA/chr_start_end.tsv

filter_file='/gpfs/commons/home/tlin/output/CARMA/chr_start_end.tsv'
input_file='/gpfs/commons/home/tlin/output/CARMA/all_variants_loc.tsv2'

# while IFS=$' ' read -r chr pos; do
#     awk -v chr="$chr" -v pos="$pos" 'BEGIN { OFS="\t"; } NR > 1 && $1 == chr && pos >= $3 && pos <= $4 { print $1, $2 }' "$filter_file"
# done < "$input_file"


while IFS=' ' read -r chr pos; do
    # Process the filter file (tab-separated) and skip the header
    match=$(awk -v chr="$chr" -v pos="$pos" '
        BEGIN { found = 0 }
        NR > 1 && $1 == chr && pos >= $3 && pos <= $4 {
            print $1, $2
            found = 1
        }
        END { if (found == 0) print "na" }
    ' "$filter_file")
    echo "$match"
done < "$input_file"