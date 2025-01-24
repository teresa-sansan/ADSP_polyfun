library(CARMA)
library(data.table)
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
chrom = args[1]
ld_id = args[2]
anno=args[3]

total_blk=3
total_blk=total_blk-1

random <- sample(letters, 2)

output_name = paste('/gpfs/commons/home/tlin/output/CARMA/',anno,'/chr', chrom, '_', ld_id,'_', sep = '')

### put the path for your annotation here
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
sumstat <- fread(file = sumstat_path , sep = "\t", header = T, check.names = F, data.table = F, stringsAsFactors = F)

ld_id <- as.numeric(ld_id)
print(sprintf('run CARMA using on chr %s, ld blk %s-%s', chrom, ld_id, ld_id + total_blk))
annot = fread(file = paste(anno_path,chrom,'.annot.gz', sep = ''),
              sep = '\t', header = T, check.names = F, data.table = F,
              stringsAsFactors = F)

z.list<-list()
ld.list<-list()
lambda.list<-list()
annot.list = list()
sumstat_sub = list() ## this is a placeholder. For annotating CARMA result.


## read data for n blk
for (i in ld_id:(ld_id+total_blk)) {
    ld_path = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/','/chr', chrom, '_', i,'.ld',sep = '' )
    ld = fread(file = ld_path, sep = " ", header = F, check.names = F, data.table = F, stringsAsFactors = F)
    snp = fread(file = paste('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/LD/LD_CARMA/geno_filt/chr', chrom, '_', i,'.bim',sep = '' ), sep = "\t", select = 2)[[1]] 
    index = i-ld_id+1
    
    sumstat_sub[[index]] <- subset(sumstat, SNP %in% snp)

    nan_count = sum(is.na(ld))
    if (nan_count > 0)
        print(sprintf('set %d nan ld score to 0', nan_count))
        ld[is.na(ld)] <- 0
    
    annot_sub <-  sumstat_sub[[index]][c('SNP','Ref','Alt')] %>%
    left_join(annot, by = "SNP") %>%
    select(SNP, everything())  

    annot_annot <- annot_sub[, -(1:7)]
    
   
    z.list[[index]] <-sumstat_sub[[index]] $Z
    ld.list[[index]] <- as.matrix(ld)
    lambda.list[[index]] <- 1   
    annot.list[[index]] = as.matrix(annot_annot)
    print(sprintf('%d snps in sumstat that overlap with LD %d',dim(sumstat_sub[[index]])[1],  i))
}

CARMA.results<-CARMA(z.list,ld.list,output_label = output_name ,lambda.list=lambda.list,w.list = annot.list, outlier.switch=T)
for (i in 1:total_blk) {
    sumstat.result = sumstat_sub[[i]] %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
    if(length(CARMA.results[[i]]$`Credible set`[[2]])!=0){
    for(l in 1:length(CARMA.results[[i]]$`Credible set`[[2]])){ sumstat.result$CS[CARMA.results[[i]]$`Credible set`[[2]][[l]]]=l} }
    output_name = paste('/gpfs/commons/home/tlin/output/CARMA/',anno,'/chr', chrom, '_', ld_id + i,'_',random, '.txt.gz', sep = '')
    fwrite(x = sumstat.result,
       file = output_name, sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
}
