
#Please include the file names and locations for each of the files below
#Default code uses the most recent (ng00067.v9) files for the 36k. To use the defaults, place each file in a folder called "data" in the same directory as this script and set working directory to source file location

cc_pheno_file = "data/ADSPCaseControlPhenotypes_DS_2022.08.18_ALL.txt"
fam_pheno_file = "data/ADSPFamilyBasedPhenotypes_DS_2022.08.18_ALL.txt"
adni_pheno_file = "data/ADNIPhenotypes_DS_2022.08.18_ALL.txt"
psp_cbd_pheno_file = "data/PSPCBDPhenotypes_DS_2022.08.18_ALL.txt"

sample_manifest_file = "data/SampleManifest_DS_2022.08.18_ALL.txt"
sample_summary_file = "data/gcad.qc.r4.wgs.allchr.36361.GATK.2022.08.15.sample.summary.ALL.txt"
ibd_file = "data/gcad.r4.wgs.36361.2022.08.15.pairwise_IBD.txt"


####################
# 0. Read data files


# Case/Control phenotypes

d.cc= read.table(cc_pheno_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")


# Family study phenotypes

d.fs= read.table(fam_pheno_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")

d.fs[d.fs=="#N/A"]=NA  # typo


# ADNI phenotypes

d.adni= read.table(adni_pheno_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")


# PSP/CBD phenotypes

d.pspcbd= read.table(psp_cbd_pheno_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")


# sample to subject mapping


sample.manifest = read.table(sample_manifest_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")

sample.manifest[sample.manifest=="#N/A"]=NA  # typo


sample.manifest=sample.manifest[sample.manifest$SAMPLE_USE %in% c("WGS"),]


# Sample Summary

sample.qc = read.table(sample_summary_file,   
                       header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")

# IBD
ibdpair = read.table(ibd_file, header=T, as.is=T,sep="\t",stringsAsFactors=FALSE,comment="",quote="")


# Define DX and Age
# AD cases use onset age and controls use last time visit ages

# Case/Control

d.cc$AGE_harmonized = NA
d.cc$DX_harmonized  = NA
d.cc$AGE_harmonized = d.cc$Age

d.cc$DX_harmonized[which(d.cc$AD==1)]=2 # AD
d.cc$DX_harmonized[which(d.cc$AD==0)]=1 # Control

d.cc$FamID = NA


# Family

d.fs$AGE_harmonized = NA
d.fs$DX_harmonized  = NA
d.fs$AGE_harmonized = d.fs$Age
d.fs$DX_harmonized[which(d.fs$AD%in%1:3)]=2 # AD
d.fs$DX_harmonized[which(d.fs$AD==0)]=1 # Control


# ADNI

d.adni$AGE_harmonized = NA
d.adni$DX_harmonized  = NA
d.adni$AGE_harmonized = NA
d.adni$DX_harmonized[which(d.adni$PrevAD == 1 | d.adni$IncAD == 1)]=2 # AD
d.adni$DX_harmonized[which(d.adni$PrevAD == 0 & d.adni$IncAD == 0 & is.na(d.adni$Age_MCI_onset) & is.na(d.adni$Age_AD_onset))]=1 # Control

d.adni$AGE_harmonized[which(d.adni $DX_harmonized==1)] = d.adni $Age_current [which(d.adni $DX_harmonized==1)] # Age for controls: current
d.adni$AGE_harmonized[which(d.adni $DX_harmonized==2)] = d.adni $Age_AD_onset[which(d.adni$DX_harmonized==2)] # Age for case: onset

d.adni$FamID = NA



####################
# 1. Identify the duplicate pairs using the threshold: PI_HAT > 0.95 from IBD file


ibdpair.95 = ibdpair[ibdpair$PI_HAT>0.95,]

ibdpair.95$SUBJID1 = sample.manifest$SUBJID[match(ibdpair.95$SampleID1,sample.manifest$SampleID)]
ibdpair.95$SUBJID2 = sample.manifest$SUBJID[match(ibdpair.95$SampleID2,sample.manifest$SampleID)]


####################
# 2. Merge the manifest file with phenotype files, drop samples from PSP/CBD studies

# subset to WGS

subj.wgs = sample.manifest$SUBJID[sample.manifest$SampleID %in% sample.qc$SampleID]

d.cc.wgs = d.cc[d.cc$SUBJID %in% subj.wgs,]
d.fs.wgs = d.fs[d.fs$SUBJID %in% subj.wgs,]
d.adni.wgs = d.adni[d.adni$SUBJID %in% subj.wgs,]
#d.pspcbd.wgs =d.pspcbd[d.pspcbd$SUBJID %in% subj.wgs,]

sample.manifest.wgs=sample.manifest[sample.manifest$SampleID %in% sample.qc$SampleID,]


####################
# 3. Drop the WGS phenotype replicates (who have phenotypes in both family and case-control studies) 
#    Rule: Keep the one in family study

n=length(intersect(d.cc.wgs$SUBJID, d.fs.wgs$SUBJID))
print(paste(c("Step 3: subjects to drop: ",n),collapse=""))

if (length(n)>0) {
  d.cc.wgs.3 = d.cc.wgs[-which(d.cc.wgs$SUBJID %in% d.fs.wgs$SUBJID),]
} else {
  d.cc.wgs.3 = d.cc.wgs
}

d.fs.wgs.3=d.fs.wgs


####################
# 4. Drop replicates (genetically identical, technique replicates) from main phenotypes

repsubj.list=unique(sample.manifest.wgs$SUBJID [sample.manifest.wgs$Technical_Replicate==1])

n=length(intersect(d.cc.wgs.3$SUBJID, repsubj.list))
print(paste(c("Step 4: subjects to drop from cc: ",n),collapse=""))


if (length(n)>0) {
  d.cc.wgs.4 = d.cc.wgs.3 [-which(d.cc.wgs.3$SUBJID %in% repsubj.list),]
} else {
  d.cc.wgs.4 = d.cc.wgs.3
}



n=length(intersect(d.fs.wgs.3$SUBJID, repsubj.list))
print(paste(c("Step 4: subjects to drop from fs: ",n),collapse=""))


if (length(n)>0) {
  d.fs.wgs.4 = d.fs.wgs.3 [-which(d.fs.wgs.3$SUBJID %in% repsubj.list),]
} else {
  d.fs.wgs.4 = d.fs.wgs.3
}



####################
# 5. combine studies into a single dataset

commonnames=names(d.cc.wgs.4);
commonnames=commonnames[commonnames%in%names(d.fs.wgs.4)];
commonnames=commonnames[commonnames%in%names(d.adni.wgs)];

d5=rbind(d.cc.wgs.4[commonnames],d.fs.wgs.4[commonnames],d.adni.wgs[commonnames])


####################
# 6.	Remove samples that failed more than one sample QC metric


# for each subject (phenotype), match with the sample with highest call rate

d5$BestSampleID=sample.manifest.wgs$SampleID[match(d5$SUBJID,sample.manifest.wgs$SUBJID)]

# update subjects with multiple samples

x=names(which(table(sample.manifest.wgs$SUBJID)>1))
x=x[x%in%d5$SUBJID]

cnt=NULL
for (subjid in x) {
  qci=sample.qc[which(sample.qc$SampleID %in% sample.manifest.wgs$SampleID[sample.manifest.wgs$SUBJID == subjid]),]
  cnt=c(cnt, nrow(qci))
  d5$BestSampleID[match(subjid,d5$SUBJID)]=qci$SampleID[which.min(qci$Missing)]
}
names(cnt)=x


####################
d7=d5

duplist.unknown = ibdpair.95[which(ibdpair.95$SUBJID1 != ibdpair.95$SUBJID2),]

duplist.unknown = duplist.unknown[ which(duplist.unknown$SUBJID1 %in% d7$SUBJID & duplist.unknown$SUBJID2 %in%  d7$SUBJID ),]


#  Define the rules to remove duplicates based on the observations in non-concordant phenotypes

resolve <- function(v1, v2) { 
  # Assume both samples are high quality
  # If APOE and AD present and inconsistent, drop
  
  v1dx=v1$DX_harmonized; if (is.na(v1dx)) {v1dx="Missing"}
  v2dx=v2$DX_harmonized; if (is.na(v2dx)) {v2dx="UNKNOWN"}
  
  v1apoe=v1$APOE_reported; if (is.na(v1apoe)) {v1apoe ="Missing"}
  v2apoe=v2$APOE_reported; if (is.na(v2apoe)) {v2apoe ="UNKNOWN"}
  
  if ((v1dx!= v2dx) || (v1apoe!= v2apoe)) { return(list("Inconsistent",NA)); }
  
  # Otherwise pick the one with more recent age
  
  # if both ages are missing, drop
  
  if (is.na(v1$AGE_harmonized) & is.na(v2$AGE_harmonized)) {return(list("BothAgesMissing",NA));}
  
  # If race/ethnicity/sex missing (unchanging data), fill with data from another sample
  
  if (is.na(v1$AGE_harmonized) & !is.na(v2$AGE_harmonized)) { select="Select1"; }
  if (!is.na(v1$AGE_harmonized) & is.na(v2$AGE_harmonized)) { select="Select2"; }
  if (!is.na(v1$AGE_harmonized) & !is.na(v2$AGE_harmonized)) {
    if (v1$AGE_harmonized>=v2$AGE_harmonized) { select="Select1"; } else { select="Select2"; }
  }
  
  # update missing data
  
  if (select == "Select1") { 
    if (is.na(v1$Race)) {v1$Race=v2$Race;}
    if (is.na(v1$Ethnicity)) {v1$Ethnicity =v2$Ethnicity;}
    if (is.na(v1$Sex)) {v1$Sex =v2$Sex;}
    return(list("Select1",v1));
  } else {
    if (is.na(v2$Race)) {v2$Race=v1$Race;}
    if (is.na(v2$Ethnicity)) {v2$Ethnicity =v1$Ethnicity;}
    if (is.na(v2$Sex)) {v2$Sex =v1$Sex;}
    return(list("Select2",v2));
  }
}

duplist.unknown$Choice=NA;  # track choices for each pair

# Go over all the duplicate pairs, merge/prune records

for (i in 1:nrow(duplist.unknown)) {
  fid1  =duplist.unknown$FID1[i]
  fid2  =duplist.unknown$FID2[i]
  subj1 =duplist.unknown$SUBJID1[i]
  subj2 =duplist.unknown$SUBJID2[i]
  rowidx1=which(d7$SUBJID==subj1)
  rowidx2=which(d7$SUBJID==subj2)
  rec1=d7[rowidx1,]
  rec2=d7[rowidx2,]
  
  if (nrow(rec1)==0) {duplist.unknown$Choice[i]="PrunedEarlier1"} else {
    if (nrow(rec2)==0) {duplist.unknown$Choice[i]="PrunedEarlier2"} else {
      
      r=resolve(rec1,rec2)
      
      duplist.unknown$Choice[i]=r[[1]]
      
      d7=d7[-c(rowidx1,rowidx2),]
      if (!(is.na(length(r[[2]]))) && !(is.na(r[[2]])))
      {
        d7=rbind(d7,r[[2]])
      }}
  }}

table(duplist.unknown$Choice)


#write.table(d7,"IntegratedPhenotype.wgs.txt",sep="\t",quote=F,row.names=F)


library(dplyr)

singles = d7$SUBJID[is.na(d7$FamID)]
families = d7[!(is.na(d7$FamID)),]


#shuffles the table to randomly select one subject (case, if available) per family
families_shuffle = families[sample(nrow(families)),]


#groups data by family, cases first, and selects the first case in each family 
one_per_fam <- 
  families_shuffle %>%
  group_by(FamID) %>%
  arrange(desc(DX_harmonized)) %>%
  filter(row_number()==1)

#combines the chosen subject per family with the individual subjects for the final subject "keep" list
one_per_fam_keep = c(singles, one_per_fam$SUBJID)

d8 = d7[(d7$SUBJID %in% one_per_fam_keep),]


write.table(d8, file="ADSP_uniqueSamples_phenotype_output.txt", row.names = FALSE, sep="\t", quote=FALSE)



