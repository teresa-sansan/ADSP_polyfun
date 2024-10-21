#cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/LD_sbayesrc
cd /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/rerun_LD_sbayesrc
rm block.eigen.bin_error2.txt
touch block.eigen.bin_error2.txt
for i in {1..2108}; do
    if [ ! -e block${i}.eigen.bin ]; then
        echo ${i} >> block.eigen.bin_error2.txt ## has 1196 in total
    fi
done
