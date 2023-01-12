1. Preprocess
ADNI (#1,565): 
- AD diagnosis is defined as AD_last_visit
- Age baseline: Age when they first entered the study
- Remove NAs 

Family_based (#4,067):
- AD diagnosis is based on the column 'AD', where:
  0, 10 -> controls (No dementia, family reported no AD)
  1-4 -> cases (Definite AD, Probable AD, Possible AD, Family Reported AD)
  5,9 -> NaN 
- Age is well-defined here, no need to specify 
- Remove NAs in APOE, Age, and Race

Case_control (# 32644)
- only need to remove the NAs

PCPCBD 
- has the AgeDeath!
- however, we cant really use this since its not AD
- PCP: Progressive supranuclear palsy
  CBD: Corticobasal degeneration


2. Merging phenotype file
Duplicates
- in Family-based and case-control
- # 112 SUBJID are duplicates, # 16 SUBJID are not just duplicates
