##############################################
# Variables: need to be fixed
ma_file="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/sbayesrc/bellenguez_hg38.cojo"     # GWAS summary in COJO format (the only input)
genotype="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/1KG/ADSP_EUR/plink/ADSP_bellenguez_snps_chr{CHR}"  # genotype prefix as LD reference (PLINK format), with {CHR} to spefify multiple
outDir="/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD_sbayesrc/LD_sbayesrc"             # Output folder that would be created automatically
threads=4                        # Number of CPU cores for eigen decomposition


#---usually don't need change bellow
genoCHR="22"                       # If more than 1 genotype file, input range (e.g. "1-22") here.
refblock=""                      # Text file to define LD blocks, by default to use our GRCH37 coordination 
tool="gctb"                      # Command line to run gctb for generating the full LD matrix


##############################################
# Code
# Step1: generate the LD block information and script
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
Rscript -e "SBayesRC::LDstep1(mafile='$ma_file', genoPrefix='$genotype', \
            outDir='$outDir', genoCHR='$genoCHR', blockRef='$refblock', log2file=TRUE)"

# Step2: generate each LD matrix for blocks
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ld.sh, $outDir/snplist/$idx.snplist
#  Ouput $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
for idx in {1..591}; do
    Rscript -e "SBayesRC::LDstep2(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
done

# Step3: eigen decomposition for each LD block
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
#  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
for idx in {1..591}; do
    Rscript -e "SBayesRC::LDstep3(outDir='$outDir', blockIndex=$idx, log2file=TRUE)"
done

# Step4: merge LD information
Rscript -e "SBayesRC::LDstep4(outDir='$outDir', log2file=TRUE)"

# Step5: clean if necessary
# Essential for analysis: $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
# Other files could be removed
# Note: before removing, check all blocks were finished.