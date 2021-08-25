#!/bin/bash
#SBATCH --job-name=bellenguez1_2.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tlin@nygenome.org
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --output=output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/%x_%j.log

summary_stats='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/bellenguez_2021/Bellenguez_2021_stage1.parquet'
bl='/gpfs/commons/groups/knowles_lab/data/ldsc/polyfun/baselineLF2.2.UKB/'
roadmap_DNase='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/DNase/roadmap_DNase_chr'
roadmap_H3K27ac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K27ac/roadmap_H3K27ac_chr'
roadmap_H3K4me1='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me1/roadmap_H3K4me1_chr'
roadmap_H3K4me3='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/roadmap_annotations/merged_annotations_ukb/H3K4me3/roadmap_H3K4me3_chr'
brain_atac='/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/annotations/brain_atac/merged_annotations_ukb/brain_atac_seq_chr'
deepsea='/gpfs/commons/home/tlin/polyfun/data/deepsea/deepsea_all_anno.'
output='output/bellenguez/bellenguez_roadmap_deepsea_brain_atac/bellenguez_roadmap_deepsea_brain_atac'

python polyfun.py \
  --compute-h2-L2 \
  --output-prefix $output \
  --sumstats $summary_stats \
  --ref-ld-chr $bl/baselineLF2.2.UKB.,$roadmap_DNase,$roadmap_H3K27ac,$roadmap_H3K4me1,$roadmap_H3K4me3,$brain_atac,$deepsea \
  --w-ld-chr $bl/weights.UKB. \
  --anno Coding_UCSC,Coding_UCSC.flanking.500,Conserved_LindbladToh,Conserved_LindbladToh.flanking.500,CTCF_Hoffman,CTCF_Hoffman.flanking.500,DGF_ENCODE,DGF_ENCODE.flanking.500,DHS_peaks_Trynka,DHS_Trynka,DHS_Trynka.flanking.500,Enhancer_Andersson,Enhancer_Andersson.flanking.500,Enhancer_Hoffman,Enhancer_Hoffman.flanking.500,FetalDHS_Trynka,FetalDHS_Trynka.flanking.500,H3K27ac_Hnisz,H3K27ac_Hnisz.flanking.500,H3K27ac_PGC2,H3K27ac_PGC2.flanking.500,H3K4me1_peaks_Trynka,H3K4me1_Trynka,H3K4me1_Trynka.flanking.500,H3K4me3_peaks_Trynka,H3K4me3_Trynka,H3K4me3_Trynka.flanking.500,H3K9ac_peaks_Trynka,H3K9ac_Trynka,H3K9ac_Trynka.flanking.500,Intron_UCSC,Intron_UCSC.flanking.500,PromoterFlanking_Hoffman,PromoterFlanking_Hoffman.flanking.500,Promoter_UCSC,Promoter_UCSC.flanking.500,Repressed_Hoffman,Repressed_Hoffman.flanking.500,SuperEnhancer_Hnisz,SuperEnhancer_Hnisz.flanking.500,TFBS_ENCODE,TFBS_ENCODE.flanking.500,Transcr_Hoffman,Transcr_Hoffman.flanking.500,TSS_Hoffman,TSS_Hoffman.flanking.500,UTR_3_UCSC,UTR_3_UCSC.flanking.500,UTR_5_UCSC,UTR_5_UCSC.flanking.500,WeakEnhancer_Hoffman,WeakEnhancer_Hoffman.flanking.500,GERP.NS,GERP.RSsup4,MAF_Adj_Predicted_Allele_Age,MAF_Adj_LLD_AFR,Recomb_Rate_10kb,Nucleotide_Diversity_10kb,Backgrd_Selection_Stat,CpG_Content_50kb,MAF_Adj_ASMC,GTEx_eQTL_MaxCPP,BLUEPRINT_H3K27acQTL_MaxCPP,BLUEPRINT_H3K4me1QTL_MaxCPP,BLUEPRINT_DNA_methylation_MaxCPP,synonymous,non_synonymous,Conserved_Vertebrate_phastCons46way,Conserved_Vertebrate_phastCons46way.flanking.500,Conserved_Mammal_phastCons46way,Conserved_Mammal_phastCons46way.flanking.500,Conserved_Primate_phastCons46way,Conserved_Primate_phastCons46way.flanking.500,BivFlnk,BivFlnk.flanking.500,Human_Promoter_Villar,Human_Promoter_Villar.flanking.500,Human_Enhancer_Villar,Human_Enhancer_Villar.flanking.500,Ancient_Sequence_Age_Human_Promoter,Ancient_Sequence_Age_Human_Promoter.flanking.500,Ancient_Sequence_Age_Human_Enhancer,Ancient_Sequence_Age_Human_Enhancer.flanking.500,Human_Enhancer_Villar_Species_Enhancer_Count,Human_Promoter_Villar_ExAC,Human_Promoter_Villar_ExAC.flanking.500,E124-DNase,E029-H3K27ac,E029-H3K4me1,E029-H3K4me3,1932,1165,1656,1168,1939,1941,microglia_atac,astrocyte_atac,neuron_atac,oligodendrocyte_atac \
  --allow-missing 




