library(ggplot2)

chirag <- read_delim("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/summary_stats/alzheimers/ADSP_reference_panel/fine_mapping/annotations_dl/aggregate_finemap_v3/combined/PRS_scores/PRS_scores_ADSP_all.tsv.gz", 
                     "\t", escape_double = FALSE, trim_ws = TRUE) %>% filter(predicted_ancestry != 'EAS')
omics_chirag  = chirag %>% filter(predicted_ancestry == 'AMR') %>% filter(cohort == 'bellenguez_omics') %>%filter(Age >=65) %>% filter(threshold <= 0.9) %>% 
  select(SampleID, threshold, SCORE_AVG, SCORESUM,Diagnosis)

omics_chirag  <- omics_chirag  %>%
  group_by(threshold) %>%
  mutate(scaled_PRS_chirag = scale(SCORE_AVG)) %>%
  mutate(scaled_SCORESUM = scale(SCORESUM))%>%
  ungroup()

AMR_teresa <- bellenguez_adsp_omics_aug %>% filter(predicted_ancestry == 'AMR') %>% select(SampleID, starts_with("PRS_"))
AMR_teresa <- AMR_teresa %>%
  rename_with(~ gsub("PRS_", "0.", .), starts_with("PRS_"))  %>%
  pivot_longer( cols = starts_with("0."),  names_to = "threshold", values_to = "PRS") %>% 
  mutate(threshold = as.numeric(threshold))%>% 
  group_by(threshold) %>%
  mutate(scaled_PRS_teresa = scale(PRS)) %>%
  ungroup()

omics_check = merge( omics_chirag , AMR_teresa,by = c("SampleID", "threshold")) %>%  mutate(threshold = as.factor(threshold))

ggplot(omics_check, aes(x = scaled_SCORESUM, y = scaled_PRS_teresa, color = threshold)) +
  geom_point() +  theme_minimal()+
  labs(
    title = "check hispanic (omics), scaled",
    x = "CHIRAG(SUM)",
    y = "TERESA",
    color = "Threshold"
  )


ggplot(omics_check, aes(x = scaled_PRS_chirag, y = scaled_PRS_teresa, color = threshold)) +
  geom_point() +  theme_minimal()+
  labs(
    title = "check hispanic, scaled (omics)",
    x = "CHIRAG(AVG)",
    y = "TERESA",
    color = "Threshold"
  )


omics_check %>%
  group_by(threshold) %>%
  filter(threshold %in% c(0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) %>%
  summarise(
    pearson_r = cor(scaled_PRS_chirag,scaled_PRS_teresa,method = "pearson"),
    .groups = "drop"
  )

omics_check %>%
  group_by(threshold) %>%
  filter(threshold %in% c(0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) %>%
  summarise(
    pearson_r = cor(scaled_SCORESUM,scaled_PRS_teresa,method = "pearson"),
    .groups = "drop"
  )

