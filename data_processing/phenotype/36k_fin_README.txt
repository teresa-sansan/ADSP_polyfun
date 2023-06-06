code is in /gpfs/commons/home/tlin/script/data_processing/phenotype/merge_pheno_36k_mergelisan.ipynb

1. Preprocess
ADNI: 
- AD diagnosis is defined as AD_last_visit
- Age: set to "Age_current" for controls, and "Age_AD_onset" for cases
- control: 958, case:608 

Family_based (#4,067):
- Age is well-defined here, no need to specify 
- AD diagnosis is based on the column 'AD', where:
  0 -> controls (No dementia) #2,493
  1-3 -> cases (Definite AD, Probable AD, Possible AD) #1,858
  others -> NaN #8,507

Case_control 
- AD diagnosis  
    #control: 23,347
    #case: 16,793
    #nan: 5,235

2. Merging phenotype file
- concatenate ADNI, Family_based, and Case control (#59,799)
- Only keep WGS (take the SUBJID found in 'gcad.qc.r4.wgs.allchr.36361.GATK.2022.08.15.sample.summary.ALL.txt') (#33,639)
- Drop tenique duplicates (#33,633)
- Duplicates in Family-based and case-control, take the ones from family_based (#33,601)

3. Adding sampleID 
    1. Merge pheotype_merge file with "SampleManifest_DS_2022.08.18_ALL.txt"
    2. For the SampleIDs that share the same SUBJID, take the one with the lowest missingess (#33,681 rows)

4. Check Kingship: Take the SampleID that overlap with the output from KING (#31,664)

5. Take one sample per family(only applies to Family_based ): (#31,517)
    - Only pick one sample from all the samples that share the same FamID. 
    - Only chose controls if there is no cases in that FamID

6. Remove rows that have NaN in Diagnosis or Age (#26,334)
    - control: 15,924
    - case: 10,410


