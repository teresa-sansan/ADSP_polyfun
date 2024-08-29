devtools::install_github("ZikunY/CARMA")
library(CARMA)

pkgs = c('data.table', 'magrittr', 'dplyr','devtools', 'R.utlis')
pkgs.na = pkgs[!pkgs %in% installed.packages()[,'Package']]
if (length(pkgs.na)>0){
  install.packages(pkgs.na)
}

library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)

setwd('/gpfs/commons/home/tlin/output/CARMA')
