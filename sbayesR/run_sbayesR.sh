#!/bin/bash

for i in 19 
do
        echo "run chr$i" 
        sbatch --export=chr=$i sbayesR.sh 

done

