library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
chrom = args[1]
ld = args[2]
anno=args[3]

maf_0.5 <- 'remove_maf_0.5/'

# if (length(args) == 3 ) {
#   maf_0.5 <- 'remove_maf_0.5/'
#   print('rerun with excluding snp with maf = 0.5')
# } else {
#   maf_0.5 <- ''
# }

if (anno == "bl") {
  anno_path <- "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations_hg38/merged_annotations_ADSP_v2/baseline_filtered/baseline_chr"
} else if (anno == "omics") {
  anno_path <- "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/teresa/anno_dec/omics_chr"
} else if (anno == "omics_dl") {
  anno_path <- "/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/teresa/anno_dec/omics_dl_chr"
} else {
  stop("Fail to find annotation")
}

#anno ='geno_filt' ## for no anno

sumstat_path = paste("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/fixed_alzheimers/processed/carma/chr",chrom,'.tsv.gz', sep = '')
ld_path = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/',maf_0.5,'/chr', chrom, '_', ld,'.ld',sep = '' )
output_name = paste('/gpfs/commons/home/tlin/output/CARMA/',anno,'/', maf_0.5,  chrom, '_', ld, '.txt.gz', sep = '')
#snp = fread(file = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt',maf_0.5,'chr', chrom, '_', ld,'.bim',sep = '' ), sep = "\t", select = 2)[[1]] 
snp = fread(file = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/chr', chrom, '_', ld,'.bim',sep = '' ), sep = "\t", select = 2)[[1]] 
print(sprintf('run CARMA using on chr %s, ld blk %s', chrom, ld))

## read data
ld = fread(file = ld_path, sep = " ", header = F, check.names = F, data.table = F, stringsAsFactors = F)
print(sprintf('%d snps in LD',dim(ld)[1]))

sumstat <- fread(file = sumstat_path , sep = "\t", header = T, check.names = F, data.table = F, stringsAsFactors = F)
sumstat <- subset(sumstat, SNP %in% snp) ## can only take the sumstat snps that's in ld blk
print(sprintf('%d snps in sumstat that overlap with LD',dim(sumstat)[1] ))

nan_count = sum(is.na(ld))
if (nan_count > 0)
  print(sprintf('set %d nan ld score to 0', nan_count))
  ld[is.na(ld)] <- 0

z.list<-list()
ld.list<-list()
lambda.list<-list()
z.list[[1]]<-sumstat$Z
ld.list[[1]]<-as.matrix(ld)
lambda.list[[1]]<-1

# ## run carma no annotation
# CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,outlier.switch=T)

# ###### Posterior inclusion probability (PIP) and credible set (CS)
# sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
# if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
#   for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
#   } }
# ###### write the GWAS summary statistics with PIP and CS
# fwrite(x = sumstat.result,
#        file = output_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")

# print(sprintf('write result in %s', output_name))


## run carma w. annotation
# annot = fread(file = paste('/gpfs/commons/home/tlin/data/ukbb_anno/baselineLF2.2.UKB.',chrom,'.annot.tsv', sep = ''),
#         sep = '\t', header = T, check.names = F, data.table = F,
#         stringsAsFactors = F)


annot = fread(file = paste(anno_path,chrom,'.annot.gz', sep = ''),
              sep = '\t', header = T, check.names = F, data.table = F,
              stringsAsFactors = F)


annot_sub <- sumstat[c('SNP','Ref','Alt')] %>%
  left_join(annot, by = "SNP") %>%
  select(SNP, everything())  

annot_annot <- annot_sub[, -(1:7)]
annot.list = list()
#annot.list[[1]] = annot_annot
annot.list[[1]] = as.matrix(annot_annot)
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list = annot.list, outlier.switch=F)

sumstat.result = sumstat %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
  for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){ sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
  } }
###### write the GWAS summary statistics with PIP and CS
fwrite(x = sumstat.result,
       file = output_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")

print(sprintf('write result in %s', output_name))