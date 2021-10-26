library(ggplot2)
library(package = "lattice")
#install.packages("sm")
library(sm)

# LOAD DATA ---------------
polypred1 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_1.tsv",header=T,fill = T)
polypred3 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_3.tsv",header=T,fill = T)
polypred5 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_5.tsv",header=T,fill = T)
polypred7 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_7.tsv",header=T,fill = T)
polypred10 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_10.tsv",header=T,fill = T)

susie_polypred1 = read.table("/gpfs/commons/home/tlin/output/prs/susie_prs_diagnosis_0219.2021_max_snp_1.tsv",header=T,fill = T)
susie_polypred3 = read.table("/gpfs/commons/home/tlin/output/prs/susie_prs_diagnosis_0219.2021_max_snp_3.tsv",header=T,fill = T)
susie_polypred5 = read.table("/gpfs/commons/home/tlin/output/prs/susie_prs_diagnosis_0219.2021_max_snp_5.tsv",header=T,fill = T)
susie_polypred7 = read.table("/gpfs/commons/home/tlin/output/prs/susie_prs_diagnosis_0219.2021_max_snp_7.tsv",header=T,fill = T)
susie_polypred10 = read.table("/gpfs/commons/home/tlin/output/prs/susie_prs_diagnosis_0219.2021_max_snp_10.tsv",header=T,fill = T)
 

polypred = rbind(polypred1, polypred3, polypred5, polypred7, polypred10)
polypred_susie = rbind(polypred1, polypred3, polypred5, polypred7, polypred10)

snp = c(1,3,5,7,10)
polypred$SNP = rep(snp,each=11903)
polypred_susie$SNP = rep(snp,each=11903)

#attach(polypred)
# --------
polypred$diagnosis_case = "Control"
polypred[polypred$diagnosis == 1,]$diagnosis_case = "Case"

polypred_susie$diagnosis_case = "Control"
polypred_susie[polypred$diagnosis == 1,]$diagnosis_case = "Case"

sort_prs <-polypred[order(PRS),]
head(sort_prs)

draw_plot <- function(PRS, percentage,main, data){
  par(mfrow=c(2,1))
  num <- round(dim(sort_prs)[1]*percentage*0.01)
  plot(head(sort_prs, num)$PRS, head(sort_prs,num)$Diagnosis,ylim = c(0,1), xlab="PRS",ylab="Diagnosis",sub = paste("top ",percentage, "%"), main=main , data = data)
  plot(tail(sort_prs, num)$PRS, tail(sort_prs,num)$Diagnosis,ylim =c(0,1),xlab="PRS",ylab="diagnosis",sub = paste("bottom ",percentage, "%"))
}

draw_plot(polypred1,10, main = "SNP = 1", data = polypred)
draw_plot(polypred1,10, main = "SNP = 10")
draw_plot(polypred, 5) ## all

control = polypred[polypred$diagnosis == 0,]$PRS
case = polypred[polypred$diagnosis == 1,]$PRS


par(mfrow=c(2,1))
hist(polypred[polypred$Diagnosis == 1,]$PRS, main = "AD patients",xlab="PRS", ylab='number', freq= FALSE)
hist(polypred[polypred$Diagnosis == 0,]$PRS, main = "AD controls",xlab="PRS", ylab='number', freq= FALSE)



par(mfrow=c(2,1))
densityplot(~ polypred$PRS, polypred[polypred$diagnosis == 1,], main = "Case")
densityplot(~ polypred$PRS, polypred[polypred$diagnosis == 0,], main = "Control")



#densityplot(~PRS, data = polypred, group = diagnosis_case, auto.key= TRUE, plot.points=FALSE)


ggplot(polypred[polypred$diagnosis == 0,]) +
  geom_histogram(aes(x = PRS,y=..density..),
                 fill = "grey", color="black", bins=30)
  #geom_line(data =  aes(x = x, y = y), color = "red"))


snp.f <- factor(polypred$SNP, levels= c(1,3,5,7,10),
                labels = c("1 snp", "3 snp", "5 snp", "7 snp", "10 snp" ))

sm.density.compare(polypred$PRS, snp.f)
colfill<-c(2:(2+length(levels(snp.f))))
legend('topright', levels(snp.f),fill=colfill)


diagnosis.f <- factor(polypred$Diagnosis, levels= c(1,0),
                labels = c("case", "control"))

sm.density.compare(polypred1$PRS, diagnosis.f)
colfill<-c(2:(2+length(levels(diagnosis.f))))
legend('topright', levels(diagnosis.f),fill=colfill)


density_plot <- function(polypred, name){
  sort = polypred[order(polypred$PRS),]$PRS
  plot(density(polypred[polypred$Diagnosis == 1,]$PRS),col = "red", main=name, xlab="PRS", xlim = c(head(sort)[6], tail(sort)[1])) 
  lines(density(polypred[polypred1$Diagnosis == 0,]$PRS), col = "blue")
}

density_plot(polypred1, "snp per locus = 1")
legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0)
density_plot(polypred10, "snp per locus = 10")
legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0)

#TEST OUT ----

par(mfrow=c(2,1))
draw_top_densityplot <- function(polypred, percentage, main){
  sort_prs <-polypred[order(polypred$PRS),]
  num <- round(dim(sort_prs)[1]*percentage*0.01)
  head <- head(sort_prs,num)
  tail <- tail(sort_prs,num)
  plot(density(head$PRS, head$Diagnosis), xlab="PRS",sub = paste("top ",percentage, "%"), main = main)
  plot(density(tail$PRS, tail$Diagnosis), xlab="PRS",sub = paste("bottom ",percentage, "%"), main =""  )
}

draw_top_densityplot(polypred1, 10, "SNP1")


par(mfrow=c(2,2))



density_plot(polypred3, "snp per locus = 3")
density_plot(polypred5, "snp per locus = 5")
density_plot(polypred7, "snp per locus = 7")
density_plot(polypred10, "snp per locus = 10")

legend("topleft", legend=c("Case", "Control"),
       col=c("red", "blue"), lty=1:1, cex=0.8,
       box.lty=0) 



# Subseting >2SD
polypred_sd = polypred[abs(scale(polypred$PRS))>2,]
polypred_susie_sd = polypred_susie[abs(scale(polypred_susie$PRS))>2,]

# check out how simlar is SuSiE vs. bellenguez w. all annotations ----
par(mfrow=c(2,2))
plot(polypred[polypred$SNP == 1,]$PRS, polypred_susie[polypred_susie$SNP == 1,]$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 1")
plot(polypred[polypred$SNP == 3,]$PRS, polypred_susie[polypred_susie$SNP == 3,]$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 3")
plot(polypred[polypred$SNP == 5,]$PRS, polypred_susie[polypred_susie$SNP == 5,]$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 5")
plot(polypred[polypred$SNP == 7,]$PRS, polypred_susie[polypred_susie$SNP == 7,]$PRS,cex = 0.1, xlab = "functional finemapping", ylab = "SuSiE", main  = "max SNP = 7")

par(mfrow=c(2,1))
plot_correlation <- function (df,main){
  plot(df[df$SNP == 1,]$PRS,  df[df$SNP == 10,]$PRS, col= ifelse(df$Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main = main)
  legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))
}
  

plot_correlation(polypred, "PolyFun")
plot_correlation(polypred_susie, "SuSiE")

## test with subset

plot_correlation(polypred_sd, "PolyFun > 2SD")
plot_correlation(polypred_susie_sd, "SuSiE >2SD ")


par(mfrow=c(2,1))
plot(polypred_sd[polypred_sd$SNP == 1,]$PRS,  polypred_susie_sd[polypred$SNP == 10,]$PRS, col= ifelse(Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main = "PolyFun")
legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))
plot(polypred_susie_sd[polypred_susie_sd$SNP == 1,]$PRS,  polypred_susie_sd[polypred_susie_sd$SNP == 10,]$PRS, col= ifelse(Diagnosis == 1, "red","black"), cex = 0.15,xlab = "max SNP = 1", ylab = "max SNP = 10", main ="SuSiE")
legend("right", legend=c("Case", "Control"),col=c("red", "black"), cex=0.8,text.font=4,inset=c(0.04,-0.0003),fill=c("red",'black'))




legend("topleft", legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8, inset=.05,fill=c("red",'black'))

# boxplot ------

par(mfrow=c(2,1))
boxplot(polypred$PRS ~ polypred$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue", main = "Bellenguez with PolyFun") 
boxplot(polypred_susie$PRS ~ polypred_susie$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue", main = "Bellenguez with SuSiE") 

plot(polypred$PRS, polypred$Diagnosis,ylim = c(0,1), xlab="PRS",ylab="diagnosis")




## multi linear regression -------

polypred_regress <- lm(Diagnosis ~ PRS + Sex + Age, data= polypred)
polypred_regress_PRS <- lm(Diagnosis ~ PRS, data= polypred)
summary(polypred_regress)
summary(polypred_regress_PRS)


polypred_susie_regress <- lm(Diagnosis ~ PRS + Sex + Age, data= polypred_susie)
anova(polypred_regress,polypred_susie)





