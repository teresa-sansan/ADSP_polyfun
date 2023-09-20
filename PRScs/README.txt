Run PRSCS
1. Check if anything in prscs.sh needs to be change. (i.e. the output dir for the log file)
2. Edit run_prscs.sh to match the file you wanna run.

Run PLINK for PRS
3. Change run_plink_36k.sh. Can either run it directly with this script or run it with run_prscs.sh, depends on if you wanna run chr chunk in parallel.
4. sum_prs.sh => will give you per PRS per individual per threshold
5. merge_prs.py => will give you a file with every PRS calculated with diff threshold, with phenotype file information.  