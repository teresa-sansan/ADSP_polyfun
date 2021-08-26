#install.packages("qqman")
#install.packages('R.utils')
#install.packages("prob")
#install.packages("heatmaply")
install.packages("pheatmap")
library("qqman")
library("R.utils")
library("prob")
library("heatmaply")
library("pheatmap")
library("ggplot2")

col_name = c("CHR","SNP","BP","A1","A2","SNPVAR","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")
par(mfrow=c(1,2)) 

## load data

#bl_max1 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl/output/finemap/max_causal_1/finemap_UKBBbaseline.extract_e-01.gz"), 
#                col.names=col_name)
#bl_max3 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl/output/finemap/max_causal_3/finemap_bl.extract_e-01.gz"), 
#                col.names=col_name)
#bl_max5 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl/output/finemap/max_causal_5/finemap_bl.extract_e-01.gz"), 
#                col.names=col_name)

#bl_susie = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl/finemap_susie/finemap_bl_susie.extract_e-01.csv.gz"), 
#                col.names=col_name)

#bl_roadmap = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/all_col/finemap/finemap_bl_roadmap.extract_e-01.gz"), 
#               col.names=col_name)
#bl_roadmap_specific_col = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/specific_col/finemap/finemap_bl_roadmap.extract_e-01.gz"), 
#                col.names=col_name)
#bl_roadmap_brainatac = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_brain_atac/finemap/finemap_bl_roadmap_extract_e-01.gz"), 
#                col.names=col_name)
#bl_roadmap_brainatac_susie = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_brain_atac/finemap_susie/finemap_bl_roadmap_brain_atac.extract_e-01.gz"), 
#                col.names=col_name)
#bl_roadmap_microglia = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_microglia/finemap/finemap_bl_roadmap_microglia.extract_e-01.gz"),
#                col.names=col_name)
# bl_roadmap_microglia_susie = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_microglia/finemap_susie/finemap_bl_roadmap_microglia.extract_e-01.gz"),
# col.names=col_name)

#bl_roadmap_susie = read.table(gzfile('/gpfs/commons/home/tlin/polyfun/output/bl_roadmap/specific_col/finemap_susie/finemap_bl_roadmap.extract_e-01.csv.gz'),
#                  col.names=col_name)


# bl_brainatac_max1 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_brain_atac/finemap/max_causal_1/finemap_UKBB_brainatac.extract_e-01.gz"), 
#                 col.names=col_name)
# 
# bl_brainatac_max3 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_brain_atac/finemap/max_causal_3/finemap_bl_brainatac_extract_e-01.gz"),
#                 col.names=col_name)
# 
# bl_brainatac_max5 = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_brain_atac/finemap/max_causal_5/finemap_bl.brain_atac.extract_e-01.gz"), 
#                 col.names=col_name)
# 
# bl_microglia = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_microglia/finemap/finemap_UKBB_brainatac.extract_e-01.gz"), 
#                 col.names=col_name)
# 
# 
# bl_brainatac_correct= read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/bl_brain_atac/finemap/correct/finemap_bl_brian_atac.extract_e-01.csv.gz"),
#                       col.names=col_name)
# bl_deepsea_brain_atac=read.table(gzfile('/gpfs/commons/home/tlin/polyfun/output/bl_deepsea_brain_atac/finemap/finemap_deepsea_brian_atac.extract_e-01.csv.gz'),
#                col.names=col_name)
# 
# bl_roadmap_deepsea = read.table(gzfile('/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_deepsea/finemap/finemap_roadmap_deepsea.extract_e-01.csv.gz'), 
#                                 col.names=col_name)
#jansen = read.table(gzfile("/gpfs/commons/home/tlin/polyfun/output/jansenetal/finemap/finemap_jansen.extract_e-01.gz"), 
#                     col.names=col_name)

bellenguez = read.table(gzfile('/gpfs/commons/home/tlin/polyfun/output/bellenguez/bellenguez/finemap/finemap_bellenguez.extract_e-01.csv.gz'))
colnames(bellenguez)<- c("CHR","SNP","BP","A1","A2","SNPVAR","MAF","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")
bellenguez = bellenguez[,-7]
bl_deepsea = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_deepsea/specific_col/finemap/finemap_bl_deepsea.extract_e-01.csv.gz',  
                      col.names=col_name)
#bl_deepsea_microglia = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_deepsea_microglia/finemap/finemap_deepsea_microglia.extract_e-01.csv.gz',
#  col.names=col_name)

#bl_pip_anno = read.table('/gpfs/commons/home/tlin/polyfun/output/bl/extract_anno/extract_pip_0.5.csv', header = TRUE,stringsAsFactors=FALSE)
#bl_microglia_pip_anno = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_microglia/extract_anno/extract_pip_0.5.csv', header = TRUE,stringsAsFactors=FALSE)

tau = read.csv('/gpfs/commons/home/tlin/polyfun/data/bl_annotation_tau.tsv', header=T, sep='\t')



# draw manhatten plot
gwas_man <- function(data, title, color=c("blue4", "firebrick1"), ylim = c(0,10),  highlight = FALSE){
    qq(data$P, main =title)
    if(highlight == FALSE){
      manhattan(data, chr = "CHR", bp = "BP", snp = "SNP", p = "P", ylim = ylim, col = color , main = title, annotatePval = 1e-8,annotateTop = TRUE)
    }
   else{
     manhattan(data, chr = "CHR", bp = "BP", snp = "SNP", p = "P", ylim = ylim, col = c("grey39","grey79") , main = title,annotatePval = 1e-8,annotateTop = FALSE, highlight = highlight)
   }
  }


 
gwas_man(brain_atac,"brain_atac")
gwas_man(microglia, "microglia")
gwas_man(data = kunkle, title = "kunkle")
gwas_man(microglia, "microglia")



qq(kunkle_gwas$P)
tail(kunkle_gwas[order(kunkle_gwas["P"]),], decreasing = FALSE)
kunkle_gwas_omit <- kunkle_gwas[!duplicated("SNP",)]


## GWAS significant

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
brain_atac_gwas = gwas_data(brain_atac)
microglia_gwas = gwas_data(microglia)
bl_kunkle_5 <- subset(baseline_kunkle, P<1e-8)
bl_kunkle_1 <- subset(baseline_kunkle1, P<1e-8)
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



##create bar plot

count_SNP <-function(df){
  return(c(dim(subset(df,P<1e-8))[1],dim(subset(df,PIP >=0.8))[1],dim(subset(df, PIP>=0.5))[1],
          dim(subset(df,PIP>=0.3))[1]))
}

create_bar_plot <- function(df,line=F){
  colnames(df) <- c("P<1e-08","PIP >= 0.8","PIP >= 0.5","PIP >= 0.3")
  only_pip <- df[,-1] #renmove P<1e-08
  barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*0.35*dim(only_pip)[1]))
  legend("topleft",rownames(df), fill=terrain.colors(dim(only_pip)[1]),bty='n')
  if(line != F)
    abline(h = c(max(test[,1]),max(test[,2]),max(test[,3])), lty=3, col = "red")
}




SNP_num <- setNames(data.frame(
  t(data.frame(count_SNP(bl_max1),count_SNP(bl_max3),count_SNP(bl_max5),
               count_SNP(bl_brainatac_max1),count_SNP(bl_brainatac_max3),count_SNP(bl_brainatac_max5),
               count_SNP(bl_microglia),count_SNP(bl_roadmap),count_SNP(bl_roadmap_specific_col),
               count_SNP(bl_roadmap_brainatac),count_SNP(bl_roadmap_microglia))),
               row.names = c('baseline(max1)', 'baseline(max3)','baseline(max5)',
                             'baseline_brain_atac(max1)', 'baseline_brain_atac(max3)','baseline_brain_atac(max5)',
                             'baseline_microglia','baseline_roadmap(all)','baselin_roadmap(specific col)',
                             'baseline_roadmap_brain_atac','baseline_roadmap_microglia')),
               c("P<1e-08","PIP>=0.8","PIP>=0.5","PIP>=0.3")
               )


SNP_main <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(bl_brainatac_max1),
                  count_SNP(bl_microglia),count_SNP(bl_roadmap_specific_col),
                  count_SNP(bl_roadmap_brainatac),count_SNP(bl_roadmap_microglia),
                  count_SNP(bl_deepsea), count_SNP(bl_deepsea_microglia))))


rownames(SNP_main) <- c('baseline','baseline_brain_atac','baseline_microglia', 'baseline_roadmap',
                        'baseline_roadmap_brain_atac','baseline_roadmap_microglia','baseline_deepsea',
                        'baseline_deepsea_microglia')
create_bar_plot(SNP_main)
create_bar_plot(SNP_main, line = T)

SNP_main_new <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(bl_brainatac_max1),
                                     count_SNP(bl_microglia),count_SNP(bl_roadmap_specific_col),
                                     count_SNP(bl_deepsea),count_SNP(bl_roadmap_deepsea))))

rownames(SNP_main_new) <- c('baseline','baseline_brain_atac','baseline_microglia', 'baseline_roadmap',
                        'baseline_deepsea','baseline_roadmap_deepsea')
create_bar_plot(SNP_main_new)

##summary_stats
SNP_summarystats <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(jansen),count_SNP(bellenguez)
)))
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

##susie vs polyfun

SNP_susie_polyfun <- data.matrix(t(
  data.frame(coun)
))





test <- SNP_num[-1]
test <- data.matrix(test[c(-2,-3,-5,-6),])

barplot(test, beside=T, width=.2, col=terrain.colors(7), ylim = c(0,68))

legend("topleft",c('baseline','baseline_brain_atac', 
      'baseline_microglia','baseline_roadmap(all)','baselin_roadmap(specific col)',
      'baseline_roadmap_brain_atac','baseline_roadmap_microglia'), fill=terrain.colors(7),bty='n')

abline(h = c(max(test[,1]),max(test[,2]),max(test[,3])), lty=3, col = "red")



## draw heat map
# 
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



# snp density between kunkle and bellenguez
common=read.csv("/gpfs/commons/home/tlin/data/common_rigorous.tsv",sep='\t',header=T)
plot(common$Z.kunkle.,common$Z.bellenguez.,xlim=c(-10,13),ylim=c(-15,22),
     xlab="kunkle",ylab="bellenguez", main="Effect Size",sub="(P < 1e-7, PIP>0.1)")
