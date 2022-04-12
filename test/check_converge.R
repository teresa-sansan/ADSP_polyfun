##check chr19 ##updataeRSID 
chr19 = read.table('bellenguez_all.19.snpvar_constrained.gz',header=T,fill = T) 
chr19_not_converge = subset(chr19, BP>43000001 & BP<48000001)
chr19_converge = subset(chr19, BP<43000001 | BP>48000001)

summary(chr19_converge$SNPVAR)
summary(chr19_not_converge$SNPVAR)


summary(chr19_converge$Z)  
summary(chr19_not_converge$Z)  


chr19_updateRSID = read.table('bellenguez_all.19.snpvar_constrained.gz',header=T,fill = T) 
chr19_fixed_0224 = read.table('/gpfs/commons/home/tlin/output/bellenguez/bellenguez_fixed_0224/bellenguez.19.snpvar_constrained.gz',header=T,fill=T)

summary(chr19_fixed_0224[1:1000,c("BP","SNPVAR","Z")])
summary(chr19_updateRSID[1:1000,c("BP","SNPVAR","Z")]) 
