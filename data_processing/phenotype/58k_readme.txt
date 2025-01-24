code is in /gpfs/commons/home/tlin/script/data_processing/phenotype/merge_pheno_58k.ipynb

1. Preprocess
ADNI(#2168): 
- AD diagnosis is defined as AD_last_visit
    - It's the same as taking either PrevAD == 1 or IncAD == 1 as case
- Age: set to "Age_current" for controls, and "Age_AD_onset" for cases
- control: 1427, case: 740 

Family_based (#42,898):
- Age is well-defined here, no need to specify 
- AD diagnosis is based on the column 'AD', where:
  0 -> controls (No dementia) # 5,022
  1-3 -> cases (Definite AD, Probable AD, Possible AD) #3,586
  others -> NaN  #34,290

Case_control (# 62,764)
- Here we use AD as AD diagnosis. (the same as using PrevAD + IncAD)
- AD diagnosis  
    #control: 33,970
    #case: 168,90
    #NaN:  11,904


2. Merge all three together (output from this step: pheno_merge.tsv)
- Only keeping the WGS. (#55,539 left)
- Checking technique replicates (None!)
- Handling duplicated SUBJID (#55,325 left)
    - Only happens between Family_based & Case_control. 199 Pairs in total.
    (1)When Diagnosis & Age are the same, use the one from Family based (remove #183)
    (2)When Diagnosis | Age were NaNs, use the another one. (remove #31 (16 *2 - 1))
        - Only 'A-LOAD-LD000442' in case_control has no NA. So only kept this one
- Check Update_ADstatus(#55,325 left)
    - if Diagnosis != 1  but Update_AD == 1, update it to be 1
    - there were 2,908 Nan and 753 controls turn into cases.
    
- Total count:
    - Diagnosis: 30,090 controls, 14,516 cases, 10,719 Nan
    - Source: 49,096 Case_control, 4061 family based, 2168 ADNI


3. Add SampleID from manifest file (output from this step: pheno_merge_w_SampleID.tsv)
    - Since SampleID and SUBJID are not one-to-one matching, we use the SampleID with the highest average.genome.coverage. 
       - This info comes from qc ("gcad.r5.wgs.58507.2024.11.03.quality.metrics.ALL.csv") which only includes SampleID
       - Manifest file is the only file that has both SampleID and SUBJID
       
       
-----------------------------------------------
To do:

# remove Kinship
# add 1KG ancestrial projection

FYI: The script for 36k:

    
    

