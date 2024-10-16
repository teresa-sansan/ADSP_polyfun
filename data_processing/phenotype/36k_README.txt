code is in /gpfs/commons/home/tlin/script/data_processing/phenotype/merge_pheno_36k.ipynb

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
- # 112 SUBJID duplicates, # 32 SUBJID of them differs
    -- they are all cases
    -- there is only one pair that the age of onset is younger in family based study (A-LOAD-LD007980, 58 yrs vs 65 yrs).
       --- there is one pair that the age of onset in case control is 18 yrs younger than the family based study (A-LOAD-LD011716, 84 yrs vs 18 yrs)
    -- 7 pairs that dont have age_baseline for case_control -> use the age of onset as age_baseline
 
 
3. Final check (total rows: #38163)
- NAN in Age_baseline: replace it with current Age (They are all cases)
- dtype: change all numeric data into int
- No SUBJID were found in "SUBJ_drop"


4. Adding sampleID
- Remove the 2 rows that appear in the SUBJ_drop file 
- SampleID and SUBJID are not 1-1 match
- ! Some SUBJID have up to 13 SampleID (which all of them are WGS)
- Extraction: 
    1. If there are both WGS and WES, take the WGS one
    2. Take either the first WGS row (#599) or the first WES row (#183).
- Some SUBJID are missing in the Manifest file. (1,917 rows)
- Ends up with 36,246 rows ('/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/ADSP_vcf/phenotype_file/release_36K/pheno_merge_sampleID.tsv')
  where there are 21,590 WGS and 14,656 WES samples
  

