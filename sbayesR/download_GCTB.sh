#!/bin/bash
#SBATCH --job-name=download_GCTB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --output=/gpfs/commons/home/tlin/data/GCTB/%x_%j.log


zenodo_get 10.5281/zenodo.3375373
zenodo_get 10.5281/zenodo.3376456
zenodo_get 10.5281/zenodo.3376628
