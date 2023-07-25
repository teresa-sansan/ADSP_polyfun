for i in {19..22}
        do
                echo "chr = $i" 
                sbatch --export=chr=$i run_polyloc.sh 
        done