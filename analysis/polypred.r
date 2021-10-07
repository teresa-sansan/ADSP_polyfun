library(ggplot2)
library(package = "lattice")
install.packages("lattice")

polypred = read.table("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/compact_filtered_vcf_16906/phenotype_data/diagnosis/prs_diagnosis_0219.2021_max_snp_10.tsv",header=T)
attach(polypred)
plot(polypred["PRS"],polypred["diagnosis"])
polypred$diagnosis_case = "Control"

polypred[polypred$diagnosis == 1,]$diagnosis_case = "Case"

sort_prs <-polypred[order(PRS),]
head(sort_prs)

draw_plot <- function(percentage){
  par(mfrow=c(2,1))
  num <- round(dim(sort_prs)[1]*percentage)
  plot(head(sort_prs, num)$PRS, head(sort_prs,num)$diagnosis,xlab="PRS",ylab="diagnosis",main = paste("top ",percentage, "%"))
  plot(tail(sort_prs, num)$PRS, tail(sort_prs,num)$diagnosis,xlab="PRS",ylab="diagnosis",main = paste("bottom ",percentage, "%"))
  
}
draw_plot(5)
draw_plot(10)

control = polypred[polypred$diagnosis == 0,]$PRS
case = polypred[polypred$diagnosis == 1,]$PRS

hist(control, main = 'controls', xlab="PRS", ylab = 'number')

ggplot(polypred[polypred$diagnosis == 0,], aes(x=case))+
         geom_histogram(aes(y=..density..))+
         geom_density(alpha=.2)

?curve

hist(polypred[polypred$diagnosis == 1,]$PRS, main = "AD patients",xlab="PRS", ylab='number', freq = FALSE)
text(x=-0.0033, y=700, labels =  round(summary(polypred[polypred$diagnosis == 0,]$PRS),4),font=1)

plot(density(case))
rug(jitter(case))


par(mfrow=c(2,1))
densityplot(~ PRS, polypred[polypred$diagnosis == 1,], main = "Case")
densityplot(~ PRS, polypred[polypred$diagnosis == 0,], main = "Control")

densityplot(~PRS|diagnosis, data = polypred,
            xlab="PRS",layout = c(1,2))


densityplot(~PRS, data = polypred, group = diagnosis_case, auto.key= TRUE, plot.points=FALSE)


ggplot(polypred[polypred$diagnosis == 0,]) +
  geom_histogram(aes(x = PRS,y=..density..),
                 fill = "grey", color="black", bins=30)
  #geom_line(data =  aes(x = x, y = y), color = "red"))


x <- seq(-0.0035,-0.0015, length.out=length(case))
with(polypred[polypred$diagnosis == 1,]$PRS,data.frame(x = x, y = dnorm(x, mean(PRS), sd(PRS))))

