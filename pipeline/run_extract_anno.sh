for i in {1..22}
do
	sbatch --export=chr=$i extract_anno.sh 
done
