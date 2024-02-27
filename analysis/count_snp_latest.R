library("R.utils")
library("prob")
library("heatmaply")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")
library("qqman")
library("dplyr")
library(UpSetR)
library(grid)
library(tidyverse)

library(mltools)
library(data.table)

install.packages('dplyr')
install.packages("tidyverse")
install.packages(c('haven','readr'))
install.packages("UpSetR")

install.packages("gridExtra")
#install.packages('tidyr', dependencies=TRUE)

col_name = c("CHR","SNP","BP","A1","A2","SNPVAR","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")
par(mfrow=c(1,1)) 

# load data ---------------------------------------------------------------
bellenguez_max_10 <- read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/test.tsv", header = T)
bellenguez_max_7 <- read.table(gzfile("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_7/agg_bellenguez.extract_1e-3.tsv.gz"), header = T, fill = T)
kunkle_max_10 <- read.table("/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10/agg_kunkle_extract_1e-3.tsv", header = T)
bellenguez_max_10_updateRSID <- read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/agg_bellenguez_extract_1e-3.tsv", header = T)
bellenguez_min_pip <- read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/agg_min_extract_1e-3.tsv", header = T)
bellenguez_max_pip <- read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_updateRSID/finemap/max_snp_10/agg_max_extract_1e-3.tsv", header = T)
bellenguez_fixed_0224_pip <- read.table('/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/finemap/max_snp_10/aggregate_extr_1e-3.tsv', header=T)
kunkle_pip = read.table("/gpfs/commons/home/tlin/output/kunkle/kunkle_fixed_0224/finemap/max_snp_10/agg_kunkle_extract_1e-3.tsv", header =T)

##wightman------------
##@ fixed_0224
wightman_susie_max5 <- read.table('/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/susie/max_snp_5/agg_extract_1e-3.tsv', header=T)
wightman_bl_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/bl/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_bl_max5 = read.table('/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/bl/max_snp_5/agg_extract_1e-3.tsv', header = T)
wightman_bl_max1 = read.table('/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/bl/max_snp_1/agg_extract_1e-3.tsv', header = T)


wightman_update_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_update_enformer_max5 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_extract_1e-3.tsv', header = T)
wightman_update_enformer_max1 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_1/agg_extract_1e-3.tsv', header = T)

## 0824 
wightman_all_max5 <- read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0824/all/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
wightman_no_ml_max5 <- read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0824/no_ml/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
wightman_only_ml_max5 <- read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0824/only_ml/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
wightman_bl_max5 <- read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0824/bl/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)

# wightman_noml_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/no_ml/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
# wightman_no_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/finemap/max_snp_10//agg_extract_1e-3.tsv', header=T)
# wightman_update_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
# wightman_enformer_max10=read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header=T)
# wightman_glasslab_max10=read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/glasslab/finemap/max_snp_10/agg_extract_1e-3.tsv', header=T)

##bellenguez
bellenguez_all_max5 <- read.table('/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/all/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
bellenguez_no_ml_max5 <- read.table('/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/no_ml/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
bellenguez_only_ml_max5 <-  read.table('/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/only_ml/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
bellenguez_bl_max5 <-  read.table('/gpfs/commons/home/tlin/output/bellenguez/new_anno_0824/bl/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)
bellenguez_susie_max5 <- read.table('/gpfs/commons/home/tlin/output/bellenguez/new_sep22/susie/finemap/max_snp_5/agg_extract_1e-3.tsv', header=T)

wightman_snp = read.table('/gpfs/commons/home/tlin/data/snp_wightman_opentarget.tsv', header = T, sep = '\t')

#tau = read.csv('/gpfs/commons/home/tlin/polyfun/data/bl_annotation_tau.tsv', header=T, sep='\t')

# visualization -----------------------------------------------------------

# *draw manhatten plot ---------
gwas_man <- function(data, title, color=c("blue4", "firebrick1"), ylim = c(5,25),  highlight = FALSE){
  qq(data$P, main =title)
  if(highlight == FALSE){
    manhattan(data, chr = "CHR", bp = "BP", snp = "SNP", p = "P", ylim = ylim, col = color , main = title, annotatePval = 1e-8,annotateTop = TRUE)
  }
  else{
    manhattan(data, chr = "CHR", bp = "BP", snp = "SNP", p = "P", ylim = ylim, col = c("grey39","grey79") , main = title,annotatePval = 1e-10,annotateTop = FALSE, highlight = highlight)
  }
}

manhattan(wightman_update_enformer_max10, annotatePval = 0.01, chr = "CHR", bp = "BP", snp = "SNP", p = "P", ylim = c(5,30), main='Wightman, 2022')


gwas_man(wightman_update_enformer_max10, "Wightman, 2022") 



# *GWAS significant ---------

gwas_data <- function(data, pvalue=FALSE, p=TRUE){
  if(pvalue == FALSE){
    pvalue=1e-8
  }
  if(p=TRUE){
    gwas_sig <- subset(data, P<pvalue)
    return(gwas_sig)
  }
  else{
    return(subset(data, p<pvalue))
  }
}

kunkle_gwas = subset(kunkle, P<1e-8)
jansen_gwas = subset(jansen, P<1e-8)
wightman_gwas = subset(wightman_update_enformer_max10, P < 1e-8)
gwas_man(wightman_gwas)

roadmap_gwas = gwas_data(roadmap)

kunkle_SNPreport = c("rs4844610","rs6733839","rs10933431","rs9271058","rs75932628","rs9473117",
                     "rs12539172","rs10808026","rs73223431","rs9331896","rs73223431","rs9331896",
                     "rs3740688","rs7933202","rs3851179","rs1128343","rs17125924","rs12881735",
                     "rs3752246","rs429358","rs6024870",
                     "rs7920721","rs138190086")


head(microglia_gwas[with(microglia_gwas,order(P)),], n=50)[1:10]
gwas_man(microglia_gwas,"microglia_gwas")

kunkle_gwas_sort <- kunkle_gwas[with(kunkle_gwas,order(P)),]

manhattan(kunkle, highlight =kunkle_SNPreport, annotatePval = 0.0005,annotateTop = FALSE, main = ("Kunkle Baseline, highlighted"))


head(baseline_kunkle[with(baseline_kunkle,order(PIP)),],n=50)[1:10]


get_report_SNP <- function(snp, data){
  init = 1
  snp_in_data = rep(NA, length(data))
  for (i in snp){
    if (isin(data$SNP,i) == TRUE){
      snp_in_data[init] = i
      init = init +1
    }
    else{
      print(i)
    }
  }
  snp_in_data <- as.character(na.omit(snp_in_data))
}


brain_atacSNP <- get_report_SNP(kunkle_SNPreport, brain_atac)
microglia_SNP <- get_report_SNP(kunkle_SNPreport, microglia)


manhattan(subset(bl_max1, PIP>0.05), p="PIP", logp = FALSE,ylim = c(0,1.1), ylab = "PIP", genomewideline = FALSE, suggestiveline = FALSE, main = "baseline PIP (max_causal = 1)")

manhattan(subset(bl_max3, PIP>0.05), p="PIP", logp = FALSE,ylim = c(0,1.1), ylab = "PIP", genomewideline = FALSE, suggestiveline = FALSE, main = "baseline PIP (max_causal = 3)")

manhattan(subset(bl_max5, PIP>0.05), p="PIP", logp = FALSE,ylim = c(0,1.1), ylab = "PIP", genomewideline = FALSE, suggestiveline = FALSE, main = "baseline PIP (max_causal = 5)")





# *bar plot (SNPs count in different PIP threshold) ------

count_SNP <-function(df){
  count = c(dim(subset(df,PIP >=0.8))[1],dim(subset(df, PIP>=0.5))[1],
            dim(subset(df,PIP>=0.3))[1])
  count_df = data.frame(
    name=c("PIP >= 0.8","PIP >= 0.5","PIP >= 0.3"),count)
  df_name = deparse(substitute(df))
  df_name = strsplit(df_name, split = '_')[[1]][-1]
  df_name =  df_name
  df_name = paste(df_name[-length(df_name)], collapse = "_")
  
  count_df$anno = df_name
  return(count_df)
}

count_SNP_new <-function(df){
  count = c(dim(subset(df,PIP >=0.8))[1],dim(subset(df, PIP < 0.8 & PIP>=0.5))[1],
            dim(subset(df,PIP < 0.5 & PIP>=0.3))[1])
  count_df = data.frame(
    name=c("PIP >= 0.8","PIP >= 0.5","PIP >= 0.3"),count)
  df_name = deparse(substitute(df))
  df_name = strsplit(df_name, split = '_')[[1]][-1]
  df_name =  df_name
  df_name = paste(df_name[-length(df_name)], collapse = "_")
  
  count_df$anno = df_name
  return(count_df)
}

wightman_0824 <- rbind(count_SNP(wightman_susie_max5),count_SNP(wightman_bl_max5), count_SNP(wightman_only_ml_max5),count_SNP(wightman_no_ml_max5),count_SNP(wightman_all_max5))

wightman_0824$anno <- factor(wightman_0824$anno, levels = c('susie','bl','no_ml','only_ml','all'))

wightman_count <- ggplot(data = wightman_0824, aes(x = anno, y = count, fill = name)) +  
  geom_bar(stat = "identity",  position = position_dodge(0.7) , alpha = 0.9, width= 0.6) +
  labs(x = "\n PIP threshold", y = "SNP count \n", title = "\n wightman_max5 \n") +  guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = count),  vjust = -0.1,  position = position_dodge(0.7), size = 3) +scale_fill_manual(values = c("skyblue","mediumturquoise", "steelblue3"))+theme_bw()

bellenguez_0824 <- rbind(count_SNP(bellenguez_susie_max5),count_SNP(bellenguez_bl_max5),count_SNP(bellenguez_only_ml_max5),count_SNP(bellenguez_no_ml_max5),count_SNP(bellenguez_all_max5))

bellenguez_0824$anno <- factor(bellenguez_0824$anno, levels = c('susie','bl','no_ml','only_ml','all'))

bellenguez_count <- ggplot(data = bellenguez_0824, aes(x = name, y = count, fill = anno)) +  
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75) +
  labs(x = "\n PIP threshold", y = "SNP count \n", title = "\n bellenguez_max5 \n") +
  geom_text(aes(label = count),  vjust = 1.5, position = position_dodge(.9), size = 3) +theme_bw() 



wightman_0824$source <- "wightman"
bellenguez_0824$source <- "bellenguez"
combined_data <- bind_rows(wightman_0824, bellenguez_0824)

# Create the facet plot
ggplot(data = combined_data, aes(x = anno, y = count, fill = name)) +  
  geom_bar(stat = "identity", position = position_dodge(0.7), alpha = 0.9, width = 0.6) +
  labs(x = NULL, y = "SNP count \n") +
  guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = count), vjust = -0.1, position = position_dodge(0.7), size = 3) +
  scale_fill_manual(values = c("skyblue", "mediumturquoise", "steelblue3")) +
  facet_wrap(~ source,scales = "free_y", ncol = 1) + theme_bw()+theme(legend.direction = "horizontal", legend.position = "bottom")      


bellenguez_0824$name <- factor(bellenguez_0824$name, levels =  c("PIP >= 0.8", "PIP >= 0.5", "PIP >= 0.3"))


## by PIP
ggplot(data = bellenguez_0824, aes(x = anno, y = count, fill = name)) +  
  geom_bar(stat = "identity",  position = position_dodge(0.7) , alpha = 0.9, width= 0.6) +
  labs(x = "\n PIP threshold", y = "SNP count \n", title = "\n bellenguez_max5 \n") +   guides(fill = guide_legend(title = NULL)) +
  geom_text(aes(label = count),  vjust = -0.1,  position = position_dodge(0.7), size = 3) +scale_fill_manual(values = c("skyblue","mediumturquoise", "steelblue3"))+theme_bw()+
  coord_cartesian(ylim=c(800,950))
#  
# bellenguez_0824_new <- rbind(count_SNP_new(bellenguez_susie_max5),count_SNP_new(bellenguez_bl_max5),count_SNP_new(bellenguez_only_ml_max5),count_SNP_new(bellenguez_no_ml_max5),count_SNP_new(bellenguez_all_max5))

# ggplot(data = bellenguez_0824_new, aes(x = anno, y = count, fill = name)) +  
#   geom_bar(stat = "identity",  position = "stack", alpha = 0.5, width = 0.6) +
#   labs(x = "\n PIP threshold", y = "SNP count \n", title = "\n bellenguez_max5 \n") +
#   geom_text(aes(label = count),  vjust = 1.5, position = "stack", size = 3) +scale_fill_manual(values = c("tan","skyblue", "steelblue3"))+theme_bw()




## SNP-COUNT-----------
wightman_susie_max5$source <- "wightman"
bellenguez_susie_max5$source <- "bellenguez"
snp_check <- bind_rows(wightman_susie_max5,bellenguez_susie_max5)

ggplot(data = snp_check, aes(x = P, color = source)) + geom_density( alpha = 0.5) +labs(x = "P value", y = "SNP count") +  theme_bw()

ggplot(data = snp_check, aes(x = PIP, color = source)) + geom_density(alpha = 0.9) +labs(x = "P value", y = "SNP count") +  theme_bw()



breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

pivot_wightman <-wightman_susie_max5 %>%
  mutate(PIP_bin = cut(PIP, breaks = breaks, labels = c("0~0.2", "0.2~0.4", "0.4~0.6", "0.6~0.8", "0.8~1"))) %>%
  group_by(PIP_bin)


# Define breaks for the bins
breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

# Create the pivot table

 wightman_susie_max5 %>%
  mutate(PIP_bin = cut(PIP, breaks = breaks, right=FALSE, labels = c("0~0.2", "0.2~0.4", "0.4~0.6", "0.6~0.8", "0.8~1"))) %>%
  group_by(PIP_bin) %>%
  summarise(count = n())

 PIP_bin = cut(wightman_susie_max5$PIP, breaks = breaks, right = FALSE, labels = c("0~0.2", "0.2~0.4", "0.4~0.6", "0.6~0.8", "0.8~1"))
 
 table(cut(wightman_susie_max5$PIP, breaks = breaks, right = FALSE))
 table(cut(bellenguez_susie_max5$PIP, breaks = breaks, right = FALSE))

