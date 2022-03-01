#!/bin/bash

for i in {19..21} 
do
        echo "run chr$i" 
        sbatch --export=chr=$i ~/script/sbayesR/sbayesR.sh 

done

