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

wightman_old = read.table("/gpfs/commons/home/tlin/output/wightman/wightman_check_1003/all_anno/finemap/max_snp_10/agg_extract_1e-3_fix_converge.tsv", header=F)
wightman_update_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_update_enformer_max5 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_5/agg_extract_1e-3.tsv', header = T)
wightman_update_enformer_max1 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_1/agg_extract_1e-3.tsv', header = T)
wightman_old = wightman_old[-17]

wightman_bl_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/fixed_0224_annotations/bl/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_noml_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/no_ml/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_no_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/all_except_enformer/finemap/max_snp_10//agg_extract_1e-3.tsv', header=T)
wightman_update_enformer_max10 = read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/update_all+enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header = T)
wightman_enformer_max10=read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/enformer/finemap/max_snp_10/agg_extract_1e-3.tsv', header=T)
wightman_glasslab_max10=read.table('/gpfs/commons/home/tlin/output/wightman/new_anno_0203/glasslab/finemap/max_snp_10/agg_extract_1e-3.tsv', header=T)

colnames(wightman_old) = colnames(wightman_update_enformer_max1)

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
  return(c(dim(subset(df,P<1e-5))[1],dim(subset(df,PIP >=0.8))[1],dim(subset(df, PIP>=0.5))[1],
          dim(subset(df,PIP>=0.3))[1]))
}
create_bar_plot <- function(df,title){
  count <- data.frame(
    name=c("PIP >= 0.8","PIP >= 0.5","PIP >= 0.3"),
    SNP=count_SNP(df)[-1]
  )
  print(count)
  plot<- barplot(height = count$SNP,names = count$name, col = 'skyblue', main = title)
  text(plot, y = count$SNP+10, labels=count$SNP, font=1, col='black')
  }

create_bar_plot(bellenguez_max_10, title="bellenguez_fixed_0224")
create_bar_plot(old_test, title="bellenguez_all_2")
create_bar_plot(bellenguez_max_10_updateRSID, title="bellenguez_updateRSID")
create_bar_plot(kunkle_max_10, title="kunkle_fixed_0224")  
create_bar_plot(wightman_update_enformer_max10, title = 'wightman, all enformer, max10')

create_bar_plot(wightman_update_enformer_max5, title = 'wightman, all enformer, max5')
create_bar_plot(wightman_update_enformer_max1, title = 'wightman, all enformer, max1')
create_bar_plot(kunkle_max_10, line=T)

create_bar_plot_legend <- function(df,remove=T,main=F){
  plot.new()
  #opar = par(oma = c(1,1,1,1), mar = c(2,2,2,12), new = TRUE)
  colnames(df) <- c("P<1e-08","PIP >= 0.8","PIP >= 0.5","PIP >= 0.3")
  only_pip <- df[,-1] #renmove P<1e-08
  if(main != F)
    barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*1.1), main=main, ylab = "SNP count")
  else
  barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*1.1), ylab = "# of SNP")
  legend("topleft",legend = c("all_enformer_max1","all_enformer_max5","all_enformer_max10"), fill=terrain.colors(dim(only_pip)[1]),cex =1.3,
         bty='n',xpd=TRUE)
  print(only_pip)
}
  
# ** different annotations ----



SNP_num <-data.matrix(data.frame(count_SNP(wightman_update_enformer_max1),count_SNP(wightman_update_enformer_max5),count_SNP(wightman_update_enformer_max10), count_SNP(wightman_old)))
#row.names(SNP_num) = c('Wightman_all(max_1)','Wightman_all(max_5)','Wightman_all(max_10)','Wightman_wo_enformer(max_10)')

create_bar_plot_legend(SNP_num)

SNP_num <- setNames(data.frame(
  t(data.frame(count_SNP(bl_max1),count_SNP(bl_max3),count_SNP(bl_max5),
               count_SNP(bl_brainatac_max1),count_SNP(bl_brainatac_max3),count_SNP(bl_brainatac_max5),
               count_SNP(bl_microglia),count_SNP(bl_roadmap),count_SNP(bl_roadmap_specific_col),
               count_SNP(bl_roadmap_brainatac),count_SNP(bl_roadmap_microglia))),
               row.names = c('baseline(max1)', 'baseline(max3)','baseline(max5)',
                             'baseline_brain_atac(max1)', 'baseline_brain_atac(max3)','baseline_brain_atac(max5)',
                             'baseline_microglia','baseline_roadmap(all)','baselin_roadmap(specific col)',
                             'baseline_roadmap_brain_atac','baseline_roadmap_microglia')),
               c("P<1e-08","PIP>=0.8","PIP>=0.5","PIP>=0.3"))


SNP_main <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(bl_brainatac_max1),
                  count_SNP(bl_microglia),count_SNP(bl_roadmap_specific_col),
                  count_SNP(bl_roadmap_brainatac),count_SNP(bl_roadmap_microglia),
                  count_SNP(bl_deepsea), count_SNP(bl_deepsea_microglia),count_SNP(kunkle_all))))



rownames(SNP_main) <- c('bl','bl_brain_atac','bl_microglia', 'bl_roadmap',
                        'bl_roadmap_brain_atac','bl_roadmap_microglia','bl_deepsea',
                        'bl_deepsea_microglia','bl_all')




SNP_num <- data.matrix(t(data.frame(count_SNP(wightman_update_enformer_max1),count_SNP(wightman_update_enformer_max5),count_SNP(wightman_update_enformer_max10), count_SNP(wightman_old))))
create_bar_plot_legend(SNP_num)


create_bar_plot(SNP_main)
create_bar_plot(SNP_main, line = T)
create_bar_plot_legend(SNP_main)


SNP_main_new <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(bl_brainatac_max1),
                                     count_SNP(bl_microglia),count_SNP(bl_roadmap_specific_col),
                                     count_SNP(bl_deepsea),count_SNP(bl_roadmap_deepsea),count_SNP(kunkle_all))))
create_bar_plot_legend(SNP_num)
rownames(SNP_main_new) <- c('baseline','baseline+brain_atac','baseline+microglia', 'baseline+roadmap',
                        'baseline+deepsea','baseline+roadmap+deepsea',"baseline+all annotations")
create_bar_plot(SNP_main_new)
create_bar_plot_legend(SNP_main_new)

# *** Replot (only PIP>0.3) ----------------
create_bar_plot_legend(SNP_main_new,main="SNP count for different annotations")

plot.new()

#opar = par(oma = c(1,1,1,1), mar = c(1,2,2,14), new = TRUE)
plt <-  barplot(SNP_main_new[,4], col=terrain.colors(dim(SNP_main_new)[1]), ylim = c(0,40), ylab = "SNP count",xaxt = "n",main ="SNP count for different annotations \n PIP ≥ 0.3")
legend(8.3,22,legend = rownames(SNP_main_new), fill=terrain.colors(dim(SNP_main_new)[1]),cex =1.1, bty='n',xpd=TRUE) 

text(plt, par("usr")[3], labels = rownames(SNP_main_new), srt = 20, adj = c(1,2), xpd = TRUE, cex=1) 


# **summary_stats ----
SNP_summarystats <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(jansen),count_SNP(bellenguez))))
rownames(SNP_summarystats) <- c("Kunkle", "Jansen", "Bellenguez")
create_bar_plot(SNP_summarystats)

bar<- barplot(SNP_summarystats[,1],  width=.2, col=terrain.colors(3), main= "p value < 1e-08")
text(bar, SNP_summarystats[,1] +85,paste("n = ",SNP_summarystats[,1]))

legend("topleft",rownames(df), fill=terrain.colors(dim(only_pip)[1]),bty='n')

SNP_selected <- data.matrix(t(
  data.frame(count_SNP(bl_max1),count_SNP(bl_max3),count_SNP(bl_max5),
             count_SNP(bl_brainatac_max1),count_SNP(bl_brainatac_max3),count_SNP(bl_brainatac_max5))))

rownames(SNP_selected) <- c('bl (max_causal_SNP=1)','bl (max_causal_SNP=3)','bl (max_causal_SNP=5)',
                            'bl+brain_atac (max_causal_SNP=1)','bl+brain_atac (max_causal_SNP=3)','bl+brain_atac (max_causal_SNP=5)')
create_bar_plot(SNP_selected)
create_bar_plot(SNP_selected,line=T)
barplot(SNP_main[,-1], beside=T, width=.2, col=terrain.colors(dim(SNP_main)[1]), ylim = c(0,max(SNP_main[,-1])))


# ** different num of max SNP per locus ----

MAX_snp <- data.matrix(t(data.frame(count_SNP(bl_all_max1),count_SNP(bl_all_max1_overlap),
                                    count_SNP(bl_all_max3),count_SNP(bl_all_max3_overlap),
                                    count_SNP(bl_all_max5),count_SNP(bl_all_max5_overlap),
                                    count_SNP(bl_all_max7), count_SNP(bl_all_max7_overlap),
                                    count_SNP(bl_all_max10),count_SNP(bl_all_max10_overlap))))

rownames(MAX_snp) <- c('max1','max1_overlap','max3','max3_overlap','max5','max5_overlap','max7','max7_overlap','max10','max10_overlap')
create_bar_plot_legend(MAX_snp,main="different max SNP per locus")  

# ** double check ---------
## load files without filtering as test
bl_max1_overlap_notuniq=read.table('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_1/anno_all.extract_e-01.csv',col.names=col_name)

subset(bl_all_max1,PIP >=0.8)
subset(bl_all_max1_overlap,PIP >=0.8)[,1:13]
subset(bl_max1_overlap_notuniq, PIP>= 0.8)[,1:13]
subset(bl_all_max3_overlap,PIP >=0.8)

## two SNPs that didnt appear in the filtered overlap file
subset(bl_max1_overlap_notuniq,SNP =='rs34971488')  ##chr 16, 
subset(bl_max1_overlap_notuniq,SNP =='rs4538760')   ##chr 6

bl_all_max1[!bl_all_max1[,2] %in% bl_max1_overlap_notuniq]


## load data with aggregration by polypred
aggregrate <- read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle/kunkle_all_2/finemap/max_snp_1/aggregrate.all.txt.gz'))

names(aggregate10) <- lapply(aggregate10[1, ], as.character)
aggregate10 <- aggregate10[-1,] 
str(aggregate10)

# bellenguez

bellenguez_all2 <- data.matrix(t(data.frame(count_SNP(bellenguez_all2_max1_overlap),count_SNP(bellenguez_all2_max3_overlap),
                                       count_SNP(bellenguez_all2_max5_overlap),count_SNP(bellenguez_all2_max7_overlap))))

rownames(bellenguez_all2) <- c('max1','max1_overlap','max3','max3_overlap','max5','max5_overlap','max7','max7_overlap','max10','max10_overlap')
create_bar_plot_legend(bellenguez_all2,main="bellenguez")  

# *susie vs polyfun -----

test <- SNP_num[-1]
test <- data.matrix(test[c(-2,-3,-5,-6),])

barplot(test, beside=T, width=.2, col=terrain.colors(7), ylim = c(0,68))

legend("topleft",c('baseline','baseline_brain_atac', 
      'baseline_microglia','baseline_roadmap(all)','baselin_roadmap(specific col)',
      'baseline_roadmap_brain_atac','baseline_roadmap_microglia'), fill=terrain.colors(7),bty='n')

abline(h = c(max(test[,1]),max(test[,2]),max(test[,3])), lty=3, col = "red")



# *heat map ------
# bl_anno_extracted <- bl_pip_extract_anno[,colSums(bl_pip_extract_anno == 0) <  round(dim(bl_pip_extract_anno)[1]*0.8)]
# bl_anno_extracted = bl_anno_extracted %>% arrange(PIP)
# row.names(bl_anno_extracted) <-  paste('Chr',bl_anno_extracted$CHR,': ',bl_anno_extracted$BP, sep = '')
# PIP = data.frame("PIP" = bl_anno_extracted$PIP)
# rownames(PIP) <- rownames(bl_anno_extracted)
# bl_anno_extracted <- bl_anno_extracted[,c(-1:-6)]
# bl_anno_extracted_mx <- as.matrix(bl_anno_extracted)
# # 
# pheatmap(bl_anno_extracted_mx,scale = 'column', main = 'Baseline PIP > 0.5',
#          cluster_cols = F, cluster_rows = F, angle_col = 315, fontsize_col = 8,
#          annotation_row = PIP)

##### heatmap ---------------------
draw_heatmap <- function(df,title){
  # omit the columns that has more than 80% of 0s
  df_extracted <- df[,colSums(df == 0) <  round(dim(df)[1]*0.8)]
  df_extracted = df_extracted %>% arrange(PIP)  # sort rows by PIP value 
  
  # extract PIP as an individual df
  row.names(df_extracted) <-  paste('Chr',df_extracted$CHR,': ',df_extracted$BP, sep = '')

  df_PIP = data.frame("PIP" = df_extracted$PIP)
  rownames(df_PIP) <- rownames(df_extracted)

  df_extracted <- df_extracted[,c(-1:-6)]  #remove cols that are not in the heatmap
  df_mx <- as.matrix(df_extracted)

  #plot
  pheatmap(df_mx, main = title,
           cluster_cols = F, scale = 'column', cluster_rows = F, angle_col = 315, fontsize_col = 8,
           annotation_row = df_PIP)
}

draw_heatmap(bl_microglia_pip_anno,'Baseline+microglia PIP > 0.5')
draw_heatmap(bl_pip_anno,"Baseline PIP > 0.5")


draw_heatmap_tau <- function(df,title){
  # omit the columns that has more than 80% of 0s
  df_extracted <- df[,colSums(df == 0) <  round(dim(df)[1]*0.8)]
  df_extracted = df_extracted %>% arrange(PIP)  # sort rows by PIP value 
  print(df_extracted)
  # extract PIP as an individual df
  row.names(df_extracted) <-  paste('Chr',df_extracted$CHR,': ',df_extracted$BP, sep = '')
  
  df_PIP = data.frame("PIP" = df_extracted$PIP)
  rownames(df_PIP) <- rownames(df_extracted)
  
  df_extracted <- df_extracted[,c(-1:-6)]  #remove cols that are not in the heatmap
  df_mx <- as.matrix(df_extracted)
  
  #plot
  pheatmap(df_mx, main = title,
           cluster_cols = F, scale = 'column', cluster_rows = F, angle_col = 315, fontsize_col = 8,
           annotation_row = df_PIP)
}

draw_heatmap_tau(bl_pip_anno,"Baseline PIP > 0.5")


# check properties ----------
# *snp density between kunkle and bellenguez ----------
common=read.csv("/gpfs/commons/home/tlin/data/common_rigorous.tsv",sep='\t',header=T)
par(mfrow=c(1,1))
plot(common$Z.kunkle.,common$Z.bellenguez.,xlim=c(-10,13),ylim=c(-15,22),
     xlab="kunkle",ylab="bellenguez", main="Effect Size",sub="(P < 1e-7, PIP>0.1)")

abline(1,1, lty=2, col="dodgerblue")


## check credible set -------
credible_1 <- subset(bellenguez_all2_1, CREDIBLE_SET > 0)
credible_3 <- subset(bellenguez_all2_3, CREDIBLE_SET > 0)
credible_5 <- subset(bellenguez_all2_5, CREDIBLE_SET > 0)
credible_7 <- subset(bellenguez_all2_7, CREDIBLE_SET > 0)
credible_10 <- subset(bellenguez_all2_10, CREDIBLE_SET > 0)

credible_1$MAX_SNP = 1
credible_3$MAX_SNP = 3
credible_5$MAX_SNP = 5
credible_7$MAX_SNP = 7
credible_10$MAX_SNP = 10


 
count_credibleset <- function(df){
  credible_vec <- vector()
  for (i in (1:10)){
    credible = dim(df[df$CREDIBLE_SET == i,])[1]
    credible_vec <- c(credible_vec, credible) 
  }
  return(credible_vec)
}

### without maxsnp = 1 -----

credible_matrix <- cbind(count_credibleset(credible_3),count_credibleset(credible_5),
                         count_credibleset(credible_7),count_credibleset(credible_10))

credible_matrix <- cbind(count_credibleset(wightman_pip), count_credibleset(wightman_update_enformer_max10), count_credibleset(wightman_update_enformer_max5), count_credibleset(wightman_update_enformer_max1))
credible_matrix



colnames(credible_matrix) <- c(3,5,7,10)



barplot(credible_matrix, xlab = "SNP counts", ylab = "max SNP per locus",
        col=brewer.pal(n = 10, name = "Paired"),main="num of SNP with credible set > 0", horiz=TRUE,xlim=c(0,3400))

#par(mar = c(0.1,0.1,0.3,0.8))
legend(3020, 4, inset=c(0,-0.0001), title="credible set",
       legend=c(1:10), fill= brewer.pal(n = 10, name = "Paired"), bty = "n", xpd=TRUE,  cex = 1.2 )

text(x = 450, y = 4.35, sprintf("total num = %s",dim(credible_10)[1]),col = "red",cex = 1.3)
text(x = 450, y = 3.1,  sprintf("total num = %s",dim(credible_7)[1]),col = "red",cex = 1.3)
text(x = 450, y = 1.9,  sprintf("total num = %s",dim(credible_5)[1]),col = "red",cex = 1.3)
text(x = 450, y = 0.76,  sprintf("total num = %s",dim(credible_3)[1]),col = "red",cex = 1.3)



### with MAX SNP = 1-----
credible_matrix <- cbind(count_credibleset(credible_1),count_credibleset(credible_3),count_credibleset(credible_5),
                         count_credibleset(credible_7),count_credibleset(credible_10))

colnames(credible_matrix) <- c(1,3,5,7,10)

barplot(credible_matrix, xlab = "SNP counts", ylab = "max SNP per locus",
        col=brewer.pal(n = 10, name = "Paired"),main="num of SNP with credible set > 0", horiz=TRUE,xlim=c(0,3400))

#par(mar = c(0.1,0.1,0.3,0.8))

legend(3020, 5.4, inset=c(0,-0.0001), title="credible set",
       legend=c(1:10), fill= brewer.pal(n = 10, name = "Paired"), bty = "n", xpd=TRUE,  cex = 1.2 )


text(x = 450, y = 5.5,  sprintf("total num = %s",dim(credible_10)[1]),col = "red",cex = 1.3)
text(x = 450, y = 4.35, sprintf("total num = %s",dim(credible_7)[1]),col = "red",cex = 1.3)
text(x = 450, y = 3.1,  sprintf("total num = %s",dim(credible_5)[1]),col = "red",cex = 1.3)
text(x = 450, y = 1.9,  sprintf("total num = %s",dim(credible_3)[1]),col = "red",cex = 1.3)
text(x = 398, y = 0.76,  sprintf("total num = %s",dim(credible_1)[1]),col = "red",cex = 1.3)







### credible set plot -------

test$POS = str_c("Chr",test$CHR,'_' ,test$start, "_",test$end)
old_test$POS = str_c("Chr",old_test$CHR,'_' ,old_test$start, "_",old_test$end)
# functions
create_PIP_subset <- function(df, thres, upperthres=TRUE){
  df <- subset(df, df$CREDIBLE_SET > 0) ## first, extract those rows that have credible set != 0
  print(names(df))
  print(str_c("Chr",df$CHR,'_' ,df$start, "_", df$end))
  df$POS = str_c("Chr",df$CHR,'_' ,df$start, "_", df$end) ## creat a new column thats easier to match in later stage
  df_count <- df[df$PIP > thres,] %>% count(POS)
  
  ## check if the boundary of PIP is a upper bound or lower bound
  if(upperthres == TRUE){
    df_count <- df[df$PIP > thres,] %>% count(POS)
  }else{
    df_count <- df[df$PIP < thres,] %>% count(POS)
  }
  
  df_count$n = as.numeric(as.character(df_count$n))  
  split = str_split_fixed(df_count$POS, "_", 3)
  df_count$CHR.BP = as.numeric(str_remove_all(str_c(split[,1],".",split[,2]),"Chr"))
  df_count = df_count[order(-df_count$CHR.BP),]
  df_count$POS <- factor(df_count$POS, levels=df_count$POS)
  return(df_count)
  
}

create_lollipop <- function(df, upperthres, lowerthres, title){
  lower <- create_PIP_subset(df, lowerthres, upperthres=FALSE)
  upper <- create_PIP_subset(df, upperthres)  
  lower$group <- paste("PIP <", lowerthres)
  upper$group <- paste("PIP >", upperthres)
  overlap <- rbind(upper,  subset(lower, lower$POS %in% upper$POS))
  
  ##drawplot
  lollipop <- ggplot(overlap, aes(x=POS, y = n), group(group)) + geom_linerange(aes(x = POS, ymin = 0, ymax = n, colour = group), position = "stack")+
    geom_point(aes(x = POS, y = n, colour = group),position = "stack")+
    theme_light() + theme_bw()+coord_flip()+
    scale_y_continuous(breaks = round(seq(min(1), max(50), by = 1),1))+
    #scale_y_continuous(breaks = round(seq(min(1), max(100), by = 5),1))+
    xlab("LD block")+ ylab("number of SNP(s)")+ggtitle(title)+
    scale_color_manual(breaks = c(paste("PIP >",upperthres),paste("PIP <", lowerthres)),
                       values=c("firebrick2", "darkblue"))+
    guides(colour = guide_legend(override.aes = list(size = 4)))
   print(lollipop)
  return(overlap)
}


PIP_0.95 <- create_lollipop(wightman_update_enformer_max10, 0.5,0.5,"Max SNP per locus = 10, kunkle")
PIP_0.5 <-create_lollipop(wightman_update_enformer_max5, 0.5,0.5,"Max SNP per locus = 5, kunkle")

#bellenguez_min_10

PIP_0.95 <- create_lollipop(bellenguez_max_pip, 0.95,0.5,"bellenguez_updateRSID, aggregate by max PIP")
PIP_0.5 <-create_lollipop(bellenguez_max_pip, 0.5,0.5,"bellenguez_updateRSID, aggregate by max PIP")

PIP_0.95 <- create_lollipop(bellenguez_min_pip, 0.95,0.5,"bellenguez_updateRSID, aggregate by min PIP")
PIP_0.5 <-create_lollipop(bellenguez_min_pip, 0.5,0.5,"bellenguez_updateRSID, aggregate by min PIP")


# newly integrate bellenguez 

PIP_0.95 <- create_lollipop(new_bellenguez_pip, 0.95,0.5,"bellenguez_fixed0224_aggregated")
PIP_0.5 <-create_lollipop(new_bellenguez_pip, 0.5,0.5,"bellenguez_fixed0224_aggregated")

PIP_0.95 <- create_lollipop(bellenguez_fixed_0224_pip, 0.95,0.5,"bellenguez_fixed_0224")

## bar plot that David asked for ---
plot_credible_bar <- function(df, title){
  df$POS = str_c("Chr",df$CHR,'_' ,df$start, "_",df$end)
  groupby <- df %>% group_by(POS, CREDIBLE_SET) 
  count_groupby <- groupby %>% count(POS) %>% filter(CREDIBLE_SET > 0)
  split = str_split_fixed(count_groupby$POS, "_", 3)
  count_groupby$CHR_BP = as.numeric(str_remove_all(str_c(split[,1],".",split[,2]),"Chr"))
  
  count_groupby = count_groupby[order(-count_groupby$CHR_BP),]
  count_groupby$CREDIBLE_SET <- factor(count_groupby$CREDIBLE_SET, levels=c(1:10))
  count_groupby$POS <- factor(count_groupby$POS, levels=unique(count_groupby$POS))
  
  ## plot them in 2 plot so that the lables won't be crushed all together. 

  first<- ggplot(data=count_groupby[(dim(count_groupby)[1] - round(dim(count_groupby)[1]/2)):dim(count_groupby)[1],] , aes(x=POS, y=n, fill=CREDIBLE_SET)) +
    geom_bar(stat="identity")+ coord_flip()+theme_light() + theme_bw() + ylab("number of SNP") + ggtitle(title)+
    scale_fill_manual(values=rep(c("#E69F00", "#56B4E9"),5)) + 
    geom_text(aes(label = n), size = 2.5, position= "stack", hjust= 0.95) 
  second <- ggplot(data=count_groupby[1:round(dim(count_groupby)[1]/2),] , aes(x=POS, y=n, fill=CREDIBLE_SET)) +
    geom_bar(stat="identity")+ coord_flip()+theme_light() + theme_bw() + ylab("number of SNP") + ggtitle(title)+
    scale_fill_manual(values=rep(c("#E69F00", "#56B4E9"),5))+
    geom_text(aes(label = n), size = 2.5, position= "stack", hjust= 0.95) 
  print(first)
  print(second)
  
  ##  plot hist of SNP count
  count_freq= as.data.frame(table(count_groupby$n))
  names(count_freq)[1] = "Credible_Set_Size"
  
  SNP_count<- ggplot(data= count_freq, aes(x = Credible_Set_Size, y = Freq)) + geom_bar(stat = "identity", fill = "skyblue") +
    theme_light() + theme_bw() + geom_text(aes(label = Freq), size = 3) +
    ggtitle(title)+xlab('The number of SNP(s) in credible set') +ylab("count")
  print(SNP_count)
  return(count_groupby)
  
}

plot_credible_bar(old_test,"Max SNP per locus = 10")
test_plot <- plot_credible_bar(old_test,"bellenguez_all_2")
test_plot <- plot_credible_bar(bellenguez_max_10,"bellenguez_fixed_0224")
plot_credible_bar(kunkle_max_10,"kunkle_fixed_0224")
plot_credible_bar(bellenguez_max_10_updateRSID, "bellenguez, updateRSID")
plot_credible_bar(bellenguez_max_pip, "bellenguez, updateRSID, aggregate by max PIP")
plot_credible_bar(bellenguez_min_pip, "bellenguez, updateRSID, aggregate by min PIP")


plot_credible_bar(bellenguez_max_10_updateRSID, "bellenguez, updateRSID")
plot_credible_bar(bellenguez_fixed_0224, "bellenguez, fixed_0224")
plot_credible_bar(new_bellenguez_pip, "bellenguez, fixed_0224, fixed convergence")


extract_SNP <- function(df){
  df$POS = str_c("Chr",df$CHR,'_' ,df$start, "_",df$end)
  df$unique = str_c(df$POS,'_' ,df$CREDIBLE_SET)
  unique_LD <-
    df%>% filter(CREDIBLE_SET>0)%>%
    group_by(unique) %>% filter(n()==1) %>% ungroup()  %>% select(-c(start,end, unique))
  print(unique_LD)
  #return(unique_LD)
}

bellenguez_updateRSID_snp = extract_SNP(bellenguez_max_10_updateRSID)
kunkle_snp =  extract_SNP(kunkle_max_10)
bellenguez_updateRSID_max_snp = extract_SNP(bellenguez_max_pip)
bellenguez_updateRSID_min_snp = extract_SNP(bellenguez_min_pip)
#bellenguez_max_10_SNP= extract_SNP(test)
#unique_LD = extract_SNP(old_test)
#kunkle_SNP = extract_SNP(kunkle_max_10)


extract_SNP(wightman_pip)

write.table(bellenguez_updateRSID_snp,"/gpfs/commons/home/tlin/data/update_RSID_bellenguez_33SNPs.tsv", row.names = FALSE, sep = '\t')

write.table(unique_LD,"/gpfs/commons/home/tlin/data/bellenguez_37SNPs.tsv", row.names = FALSE, sep = '\t')

par(mfrow=c(1,1))
plot(x = bellenguez_updateRSID_snp$CHR, y=bellenguez_updateRSID_snp$PIP, xlab = "CHR", ylab = "PIP", main= 'the only SNP in a credible set (bellenguez)',axes = FALSE, xlim =c(1,22))
box(bty = "l")
axis(1 , at=c(1:22))
axis(2, at =c(0,0.25,0.5,0.75,1))
lowpip <- bellenguez_updateRSID_snp %>%
  filter(PIP < 0.5) 

text(as.numeric(unlist(lowpip[,'CHR'])), as.numeric(unlist(lowpip[,'PIP'])),labels=(as.character(lowpip$SNP)),  
     cex=1, pos=3,col="blue") 
interested_SNP <- bellenguez_updateRSID_snp %>%
  filter(PIP > 0.5) 

interested_SNP_threshold <- bellenguez_max_10_updateRSID %>%
  filter(CREDIBLE_SET != 0) %>%
  filter(PIP > 0.5)

interested_SNP_kunkle <- kunkle_pip %>%
  filter(PIP > 0.5)

interested_SNP_wightman <- wightman_pip %>%
  filter(PIP > 0.95)


interested_SNP_bellenguez_new <- new_bellenguez_pip %>%
  filter(PIP > 0.5)


interested_SNP_bellenguez_notfixconvergence <- bellenguez_fixed_0224_pip %>%
  filter(CREDIBLE_SET != 0) %>%
  filter(PIP > 0.5)


write.table(interested_SNP,"/gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP.tsv", row.names = FALSE, sep = '\t',quote=F)

write.table(bellenguez_updateRSID_min_snp[bellenguez_updateRSID_min_snp$PIP> 0.5,], "/gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP_minPIP.tsv", row.names = FALSE, sep = '\t',quote=F)
write.table(bellenguez_updateRSID_max_snp[bellenguez_updateRSID_max_snp$PIP> 0.5,], "/gpfs/commons/home/tlin/data/bellenguez_updateRSID_interested_SNP_maxPIP.tsv", row.names = FALSE, sep = '\t',quote=F)

bellenguez_updateRSID_snp[bellenguez_updateRSID_snp$PIP> 0.95,]$SNP
bellenguez_updateRSID_min_snp[bellenguez_updateRSID_min_snp$PIP> 0.95,]$SNP
bellenguez_updateRSID_max_snp[bellenguez_updateRSID_max_snp$PIP> 0.95,]$SNP

## another way of thinking?
kunkle_max_10$pos = str_c("Chr",kunkle_max_10$CHR,'_' ,kunkle_max_10$start, "_",kunkle_max_10$end, '_', kunkle_max_10$CREDIBLE_SET)
otherway_kunkle = kunkle_max_10 %>% filter(PIP >=0.5 ) %>% group_by(pos) %>% ungroup()
check_kunkle = kunkle_max_10 %>% group_by(pos) %>% filter(n()==1) %>% ungroup()  %>% filter(PIP >0 )
 
bellenguez_max_10_updateRSID$pos = str_c("Chr",bellenguez_max_10_updateRSID$CHR,'_' ,bellenguez_max_10_updateRSID$start, "_",bellenguez_max_10_updateRSID$end, '_', bellenguez_max_10_updateRSID$CREDIBLE_SET)
dim(bellenguez_max_10_updateRSID %>% filter(PIP >0))  ##37062
check_bellenguez_snp <-
  bellenguez_max_10_updateRSID %>% group_by(pos) %>% filter(n()==1) %>% ungroup()  %>% filter(PIP >0 ) ##106


bellenguez_max_10_updateRSID %>% group_by(pos) %>% filter(n()==1) %>% ungroup()  %>% filter(PIP >0.5 ) %>% filter(CREDIBLE_SET ==0)

check_ld = bellenguez_max_10_updateRSID[bellenguez_max_10_updateRSID$pos == 'Chr17_43000001_46000001_1',]

write.table(check_ld,"/gpfs/commons/home/tlin/data/check_bellenguez_ld.tsv", row.names = FALSE, sep = '\t',quote=F)
write.table(check_ld$SNP,"/gpfs/commons/home/tlin/data/check_bellenguez_ld_snp_only.tsv", row.names = FALSE, sep = '\t',quote=F)




## COUNT SNP------
library('dplyr')

wightman_update_enformer_max1 %>%
  dplyr::filter(PIP > 0.9)

wightman_update_enformer_max1[wightman_update_enformer_max1$PIP > 0.9,]
upset_SNP <- function(df, filter){
  PIP8 <- df[df$PIP > 0.8,]
  PIP8$PIP_FILT = 0.8
  PIP5 <- df[df$PIP > 0.5,]
  PIP5$PIP_FILT = 0.5
  PIP3 <- df[df$PIP > 0.3,]
  PIP3$PIP_FILT = 0.3
  return(c(PIP8, PIP5, PIP3))
  
}



SNP_filter <- function(df, pip_filter){
  df<- df[df$PIP > pip_filter,]
  df$PIP_FILT = pip_filter
  return(df)
}


Upset_SNP <- function(df,anno_name,filt){
  ## create the dataset
  PIP <- SNP_filter(df, filt)
  new_df = PIP
  new_df["ANNO"] = anno_name
  new_df = new_df[,c("CHR","SNP","BP","SNPVAR","P","PIP","BETA_MEAN","BETA_SD","PIP_FILT","ANNO")]
  return(new_df)
}

Upset_SNP_plot <- function(filt){
  all <- Upset_SNP(wightman_update_enformer_max10,"all_anno",filt)
  bl  <- Upset_SNP(wightman_bl_max10,'baseline',filt) 
  all_except_enformer <- Upset_SNP(wightman_no_enformer_max10,'no_enformer',filt)
  no_ml <- Upset_SNP(wightman_no_ml_10,"no_ml",filt)
  enformer <- Upset_SNP(wightman_enformer_max10, "enformer",filt)
  glasslab <- Upset_SNP(wightman_glasslab_max10, "glasslab",filt)
  try = rbind(all, bl, all_except_enformer, no_ml, enformer, glasslab)
 
  try['bl'] = 1 ## bc every anno has bl
  try  <- mutate(try , roadmap_deepsea = ifelse(ANNO=="all_anno" | ANNO == 'no_enformer',1,0))
  try  <- mutate(try , enformer = ifelse(ANNO=="all_anno"|ANNO == 'enformer',1,0))
  try  <- mutate(try ,  glass_lab = ifelse(ANNO=="all_anno" | ANNO== 'no_enformer' | ANNO== 'no_ml'|ANNO =='glasslab' ,1,0))
  try  <- mutate(try , glasslab_enformer = ifelse(ANNO=="all_anno" |ANNO =='glasslab' |ANNO == 'enformer',1,0))

  print(upset(try, sets = c("bl", "roadmap_deepsea", "glass_lab", "enformer","glasslab_enformer"), sets.bar.color = "#56B4E9",order.by = "freq",keep.order = TRUE,set_size.show = FALSE))
  return(try)
}

pip_0.8 <- Upset_SNP_plot(0.8)
grid.text("PIP>0.8",x = 0.65, y=0.95, gp=gpar(fontsize=15))

pip_0.5 <- Upset_SNP_plot(0.5)
grid.text("PIP>0.5",x = 0.65, y=0.95, gp=gpar(fontsize=15))

pip_0.3 <- Upset_SNP_plot(0.3)
grid.text("PIP>0.3",x = 0.65, y=0.95, gp=gpar(fontsize=15))

df = rbind(pip_0.8, pip_0.5, pip_0.3)
df$PIP_FILT <- as.factor(df$PIP_FILT)
barplot(rbind(pip_0.8, pip_0.5, pip_0.3), beside=T, 
        col=c("aquamarine3","coral"), 
        names.arg=LETTERS[1:2])


ggplot(rbind(pip_0.8, pip_0.5, pip_0.3), aes(factor(ANNO, level=c('all_anno', 'enformer', 'no_ml','glasslab','no_enformer','baseline')), fill=factor(PIP_FILT))) + 
  geom_bar(position="dodge") + theme_bw()+ scale_fill_discrete(name = "PIP FILTER")+xlab('') + ylab('SNP count')+ coord_flip()+theme(legend.position="bottom")

