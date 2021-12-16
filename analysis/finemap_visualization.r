library("R.utils")
library("prob")
library("heatmaply")
library("pheatmap")
library("ggplot2")
library("RColorBrewer")
library("dplyr")
col_name = c("CHR","SNP","BP","A1","A2","SNPVAR","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")
par(mfrow=c(1,2)) 

# load data ---------------------------------------------------------------


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

# kunkle_all=read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap/all_anno_extract_e-01.csv.gz'),
                      # col.names=col_name)

# bl_all_max1 = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_1/all_anno_extract_e-01.csv.gz'), col.names=col_name)
# bl_all_max3 = read.table('/gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_3/anno_all.extract_e-01.csv', col.names=col_name)
# bl_all_max5 = read.table('/gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_5/anno_all.extract_e-01.csv', col.names=col_name)
# bl_all_max7 = read.table('/gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_7/anno_all.extract_e-01.csv', col.names=col_name)
# bl_all_max10 = read.table('/gpfs/commons/home/tlin/output/kunkle_all/finemap/finemap_max_snp_10/anno_all.extract_e-01.csv', col.names=col_name)
# 
# 
# bl_all_max1_overlap = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_1/all_anno_extract_uniq_e-01.csv.gz'), col.names=col_name)
# bl_all_max3_overlap = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_3/all_anno_extract_uniq_e-01.csv.gz'), col.names=col_name)
# bl_all_max5_overlap = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_5/all_anno_extract_uniq_e-01.csv.gz'), col.names=col_name)
# bl_all_max7_overlap = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_7/all_anno_extract_uniq_e-01.csv.gz'), col.names=col_name)
# bl_all_max10_overlap = read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all/finemap_overlap/finemap_max_snp_10/all_anno_extract_uniq_e-01.csv.gz'), col.names=col_name)

col_name_bellenguez =  c("CHR","SNP","BP","A1","A2","SNPVAR","MAF","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")

# bellenguez = read.table(gzfile('/gpfs/commons/home/tlin/polyfun/output/bellenguez/bellenguez/finemap/finemap_bellenguez.extract_e-01.csv.gz'))
colnames(bellenguez)<- c("CHR","SNP","BP","A1","A2","SNPVAR","MAF","N","Z","P","PIP","BETA_MEAN","BETA_SD","DISTANCE_FROM_CENTER","CREDIBLE_SET")
bellenguez = bellenguez[,-7]
bl_deepsea = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_deepsea/specific_col/finemap/finemap_bl_deepsea.extract_e-01.csv.gz',  
                      col.names=col_name)
#bl_deepsea_microglia = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_deepsea_microglia/finemap/finemap_deepsea_microglia.extract_e-01.csv.gz',
#  col.names=col_name)

#bl_pip_anno = read.table('/gpfs/commons/home/tlin/polyfun/output/bl/extract_anno/extract_pip_0.5.csv', header = TRUE,stringsAsFactors=FALSE)
#bl_microglia_pip_anno = read.table('/gpfs/commons/home/tlin/polyfun/output/bl_roadmap_microglia/extract_anno/extract_pip_0.5.csv', header = TRUE,stringsAsFactors=FALSE)


bellenguez_susie1 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/max_snp_1/finemap_bellenguez_susie.extract_1e-3.csv.gz",header = TRUE)
bellenguez_susie3 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_bl/finemap_susie/max_snp_3/finemap_bellenguez_susie.extract_1e-3.csv.gz", header = TRUE)

bellenguez_all2_1 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_1/finemap_bellenguez_all_2.extract_1e-3.csv.gz",header = TRUE )
bellenguez_all2_3 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_3/finemap_bellenguez_all_2.extract_1e-3.csv.gz", header =TRUE)
bellenguez_all2_5 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_5/finemap_bellenguez_all_2.extract_1e-3.csv.gz",header = TRUE )
bellenguez_all2_7 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_7/finemap_bellenguez_all_2.extract_1e-3.csv.gz", header =TRUE)
bellenguez_all2_10 = read.table("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/finemap_bellenguez_all_2.extract_1e-3.csv.gz",header = TRUE )


tau = read.csv('/gpfs/commons/home/tlin/polyfun/data/bl_annotation_tau.tsv', header=T, sep='\t')


aggregrate10 = read.table(gzfile("/gpfs/commons/home/tlin/output/bellenguez/bellenguez_all_2/finemap_snpvar_constrained/max_snp_10/finemap_bellenguez_all_2.extract_1e-3.csv.gz"), header = T)



# visualization -----------------------------------------------------------

# *draw manhatten plot ---------
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



# *bar plot (SNPs count in different categories) ------

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


create_bar_plot_legend(SNP_main_new)

create_bar_plot_legend <- function(df,remove=T,main=F){
  plot.new()
  #opar = par(oma = c(1,1,1,1), mar = c(2,2,2,12), new = TRUE)
  colnames(df) <- c("P<1e-08","PIP >= 0.8","PIP >= 0.5","PIP >= 0.3")
  only_pip <- df[,-1] #renmove P<1e-08
  if(main != F)
    barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*1.1), main=main, ylab = "SNP count")
  else
  barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*1.1), ylab = "SNP count")
  #barplot(only_pip, beside=T, width=.2, col=terrain.colors(dim(only_pip)[1]), ylim = c(0,max(only_pip)*1.1))
  
  legend(0.3,35,legend = rownames(df), fill=terrain.colors(dim(only_pip)[1]),cex =1.3,
         bty='n',xpd=TRUE)
  
  par(opar)
  print(only_pip)
    
}
  



# ** different annotations ----

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
                  count_SNP(bl_deepsea), count_SNP(bl_deepsea_microglia),count_SNP(kunkle_all))))


rownames(SNP_main) <- c('bl','bl_brain_atac','bl_microglia', 'bl_roadmap',
                        'bl_roadmap_brain_atac','bl_roadmap_microglia','bl_deepsea',
                        'bl_deepsea_microglia','bl_all')
create_bar_plot(SNP_main)
create_bar_plot(SNP_main, line = T)
create_bar_plot_legend(SNP_main)


SNP_main_new <- data.matrix(t(data.frame(count_SNP(bl_max1),count_SNP(bl_brainatac_max1),
                                     count_SNP(bl_microglia),count_SNP(bl_roadmap_specific_col),
                                     count_SNP(bl_deepsea),count_SNP(bl_roadmap_deepsea),count_SNP(kunkle_all))))
create_bar_plot_legend(SNP_main_new)
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
aggregrate <- read.table(gzfile('/gpfs/commons/home/tlin/output/kunkle_all_2/finemap/max_snp_1/aggregrate.all.txt.gz'))



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

# functions
create_PIP_subset <- function(df, thres, upperthres=TRUE){
  df <- subset(df, df$CREDIBLE_SET > 0) ## first, extract those rows that have credibleset != 0
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
    xlab("LD block")+ ylab("number of SNP(s)")+ggtitle(title)+
    scale_color_manual(breaks = c(paste("PIP >",upperthres),paste("PIP <", lowerthres)),
                       values=c("firebrick2", "darkblue"))+
    guides(colour = guide_legend(override.aes = list(size = 4)))
   print(lollipop)
  return(overlap)
}


PIP_0.95 <- create_lollipop(aggregrate10, 0.95,0.5,"Max SNP per locus = 10")
PIP_0.5 <-create_lollipop(aggregrate10, 0.5,0.5,"Max SNP per locus = 10")



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
  
}

plot_credible_bar(aggregrate10,"Max SNP per locus = 10")
