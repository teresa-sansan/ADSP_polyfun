path=/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_sbayesrc/ADSP_hg38/LD_sbayesrc
cd $path
touch ../rerun_LD_sbayesrc/rerun_blk.txt

for i in {1..2096}
do
if [ ! -e b${i}.ldm.full.bin ];then
    echo $i | tee -a ../rerun_LD_sbayesrc/rerun_blk.txt
fi
done
