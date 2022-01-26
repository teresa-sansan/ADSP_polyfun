# ** tutorial ------
height_change_target <- read.delim("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_target_file/height_prs")
phenotype <- read.table("/gpfs/commons/home/tlin/data/plink_tutorial/EUR.height", header = T)
height_pthres <- c(0.5,0.1,0.05,0.01,0.001,0.005,0.3,0.2)
pcs <- read.table("/gpfs/commons/home/tlin/data/plink_tutorial/EUR.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 

covariate <- read.table("/gpfs/commons/home/tlin/data/plink_tutorial/EUR.cov", header=T)
# Now merge the files
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c('FID','IID'))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as height is quantitative)
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# And the R2 of the null model is 
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL
for(i in height_pthres){
  # Go through each p-value threshold
  # Merge the prs with the phenotype matrix
  # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
  # relevant columns
  pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
  # Now perform a linear regression on Height with PRS and the covariates
  # ignoring the FID and IID from our model
  model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
  # model R2 is obtained as 
  model.r2 <- summary(model)$r.squared
  # R2 of PRS is simply calculated as the model R2 minus the null R2
  prs.r2 <- model.r2-null.r2
  # We can also obtain the coeffcient and p-value of association of PRS as follow
  prs.coef <- summary(model)$coeff["SCORE",]
  prs.beta <- as.numeric(prs.coef[1])
  prs.se <- as.numeric(prs.coef[2])
  prs.p <- as.numeric(prs.coef[4])
  # We can then store the results
  prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]




# ** change base data ------
library(ggplot2)
height_tutorial <- read.delim("/gpfs/commons/home/tlin/data/plink_tutorial/height_prs.tsv")
height_change_base <- read.delim("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_base_file/height_prs.tsv")
height_change_target <- read.delim("/gpfs/commons/home/tlin/data/plink_tutorial/height_test/change_target_file/height_prs")

height_tutorial["source"] = "tutorial"
height_change_base["source"] = "change_base"
height_change_target["source"] = "change_target"
df <- height_change_target[,c("source","prs_001","prs_1","prs_5")]

colnames(df)  <- c("source","PRS_001","PRS_1","PRS_5")

prs_height <- rbind(height_change_base[,c("source","PRS_001","PRS_1","PRS_5")],
                    df,
                   height_tutorial[,c("source","PRS_001","PRS_1","PRS_5")]) 


ggplot(prs_height, aes(x=PRS_001)) + geom_density()+
  facet_grid(source ~. ,scales = "free")+
  ggtitle("Height (sanity check), pT = 0.001")

ggplot(prs_height, aes(x=PRS_5)) + geom_density()+
  facet_grid(source ~. ,scales = "free")+
  ggtitle("Height (sanity check), pT = 0.5")


par(mfrow=c(2,2))
plot(density(height_change_base$PRS_001), main = "0.001")

ggplot(height_change_base, aes(x=PRS_001)) + geom_density()+
  geom_vline(aes(xintercept=mean(PRS_001, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

ggplot(height_change_base, aes(x=PRS_001)) + geom_density()




