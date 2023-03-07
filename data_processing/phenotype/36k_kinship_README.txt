This file explain how I handle the duplicated SUBJID when assiging sampleID. (multiple sampleID to 1 SUBJID situation)
Code is in /gpfs/commons/home/tlin/script/data_processing/phenotype/merge_pheno_36k_kinship.ipynb
All the files are in /gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K

Directly starts with
1) The merged phenotype file: phenotype_merged.tsv (This is the integration of three phenotype files, before adding SampleID.)
2) Non-related individual's SampleID file After running King (Thanks Anjali)): unrelated_samples.txt
3) Manifest (The file that consists both SUBJID and sampleID): SampleManifest_DS_2022.08.18_ALL.txt


Process(For Manifest file):
1. Remove the 2 rows that appear in the SUBJ_drop file 

2. Take the intersection with unrelated_samples.txt (If the SampleID is found in unrelated_samples.txt, assign that SampleID to the SUBJID.)
    3,987 out of the 4,306 duplicated SUBJID were found in unrelated_samples.txt 

3. For the duplicated SUBJIDs that were not in unrelated_samples.txt (#342 SUBJID; #735 rows):
    - If there are both WGS and WES, take the WGS one
    - If there are still duplicates, take either the first WGS row or the first WES row
    - end up with 154 WGS and 129 WES 
    
4. Merge the merged_manifest file with merged_phenotype file
- Ends up with 36,210 rows (output file: pheno_merge_sampleID_king.tsv)
  where there are 21,608 WGS and 14,602 WES samples
  

