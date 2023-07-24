for i in {1..21}
        do
                echo "chr = $i" 
                sbatch --export=chr=$i run_polyloc.sh 
        done