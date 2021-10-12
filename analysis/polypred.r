library(ggplot2)
library(package = "lattice")
#install.packages("sm")
library(sm)

# LOAD DATA ---------------
polypred1 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_1.tsv",header=T)
polypred3 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_3.tsv",header=T)
polypred5 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_5.tsv",header=T)
polypred7 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_7.tsv",header=T)
polypred10 = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_10.tsv",header=T)

polypred1$SNP = 1
polypred3$SNP = 3
polypred5$SNP = 5
polypred7$SNP = 7
polypred10$SNP = 10

polypred = rbind(polypred1, polypred3, polypred5, polypred7, polypred10)
attach(polypred)
# --------
polypred$diagnosis_case = "Control"
polypred[polypred$diagnosis == 1,]$diagnosis_case = "Case"

sort_prs <-polypred[order(PRS),]
head(sort_prs)

draw_plot <- function(PRS, percentage,main){
  par(mfrow=c(2,1))
  num <- round(dim(sort_prs)[1]*percentage*0.01)
  plot(head(sort_prs, num)$PRS, head(sort_prs,num)$diagnosis,ylim = c(0,1), xlab="PRS",ylab="diagnosis",sub = paste("top ",percentage, "%"), main=main )
  plot(tail(sort_prs, num)$PRS, tail(sort_prs,num)$diagnosis,ylim =c(0,1),xlab="PRS",ylab="diagnosis",sub = paste("bottom ",percentage, "%"))
}

draw_plot(polypred1,10, main = "SNP = 1")
draw_plot(polypred1,10, main = "SNP = 10")
draw_plot(polypred, 5) ## all

control = polypred[polypred$diagnosis == 0,]$PRS
case = polypred[polypred$diagnosis == 1,]$PRS



hist(polypred[polypred$Diagnosis == 1,]$PRS, main = "AD patients",xlab="PRS", ylab='number', freq= FALSE)
hist(polypred[polypred$Diagnosis == 0,]$PRS, main = "AD controls",xlab="PRS", ylab='number', freq= FALSE)




par(mfrow=c(2,1))
densityplot(~ PRS, polypred[polypred$diagnosis == 1,], main = "Case")
densityplot(~ PRS, polypred[polypred$diagnosis == 0,], main = "Control")



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


# boxplot ------
par(mfrow=c(1,1))
boxplot(polypred$PRS ~ polypred$SNP, horizontal=TRUE, xlab="PRS", ylab = "num of SNP per locus", col="dodgerblue") 






