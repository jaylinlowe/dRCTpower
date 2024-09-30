#' This is mostly a copy of 00-data-setup.Rmd from https://github.com/manncz/edm-rct-tutorial/tree/main/scripts
#' Info on preprocessing of the dataset loaded here can be found on https://github.com/manncz/aeis-aux-rct/

library(optmatch)
library(dplyr)
library(tidyr)
library(kableExtra)
library(randomForest)



load("MS_data_public.Rdata")
var.names.clean <- read.csv("var_names.csv")


## Clean some variable names for interpretability
covs_ms <- covs_ms %>%
  rename(all_stud_n = CPETALLC_67, grade8_n = CPETG08C_67, stud_teach_rat = CPSTKIDR_67,
         all_exp = CPFPAALLT_67, inst_exp = CPFEAINSP_67, lead_exp = CPFEAADSP_67, supp_exp = CPFEASUPP_67,
         ed_aide = CPSETOFC_67, teach_salary = CPSTTOSA_67,
         teach_expr = CPSTEXPA_67, perc_teach_female = CPSTFEFP_67,
         perc_teach_white = CPSTWHFP_67, perc_teach_black = CPSTBLFP_67, perc_teach_hisp = CPSTHIFP_67,
         perc_stud_black = CPETBLAP_67, perc_stud_hisp = CPETHISP_67,
         perc_stud_api = CPETPACP_67, perc_stud_white = CPETWHIP_67,
         perc_stud_alp = CPETSPEP_67, perc_stud_bil = CPETBILP_67, perc_stud_tag = CPETGIFP_67) %>%
  select(CAMPUS, GRDSPAN, starts_with("pre"),all_of(var.names.clean$var_clean), everything())

covs_ms_noprep <- covs_ms_noprep %>%
  rename(all_stud_n = CPETALLC_67, grade8_n = CPETG08C_67, stud_teach_rat = CPSTKIDR_67,
         all_exp = CPFPAALLT_67, inst_exp = CPFEAINSP_67, lead_exp = CPFEAADSP_67, supp_exp = CPFEASUPP_67,
         ed_aide = CPSETOFC_67, teach_salary = CPSTTOSA_67,
         teach_expr = CPSTEXPA_67, perc_teach_female = CPSTFEFP_67,
         perc_teach_white = CPSTWHFP_67, perc_teach_black = CPSTBLFP_67, perc_teach_hisp = CPSTHIFP_67,
         perc_stud_black = CPETBLAP_67, perc_stud_hisp = CPETHISP_67,
         perc_stud_api = CPETPACP_67, perc_stud_white = CPETWHIP_67,
         perc_stud_alp = CPETSPEP_67, perc_stud_bil = CPETBILP_67, perc_stud_tag = CPETGIFP_67) %>%
  select(CAMPUS, GRDSPAN, starts_with("pre"),all_of(var.names.clean$var_clean), everything())

## Form Fake RCT

set.seed(25)

poss_rct <- covs_ms_noprep %>%
  filter(!is.na(grade8_n)) %>%
  select(CAMPUS) %>%
  left_join(covs_ms)

# create a random grouping to be able to do bipartite matching
poss_rct$Z <- sample(0:1, nrow(poss_rct), replace = T)

# calculate distances between schools
dist_mat <- match_on(Z ~ premA + premB + premF + all_stud_n + stud_teach_rat + all_exp +
                       perc_teach_white + perc_stud_hisp +
                       perc_stud_alp + perc_stud_bil, data = poss_rct, caliper = 1)

# exclude distances outside of the caliper
initialmatch <- fullmatch(dist_mat, data = poss_rct)
excl <- which(is.na(initialmatch))

# and recalculate distances
dist.update <- match_on(Z ~ premA + premB + premF + all_stud_n + stud_teach_rat + all_exp +
                          perc_teach_white + perc_stud_hisp +
                          perc_stud_alp + perc_stud_bil, dat = poss_rct[-excl,])

# final optimal pair matching
fakematch <- pairmatch(dist.update,  data = poss_rct)
summary(fakematch)


# select 25 random pairs from those matched to make up our paired trial.


set.seed(125)

match.ids <- unique(fakematch[which(!is.na(fakematch))])
rct.ids <- sample(match.ids, 25)


# Randomly assign one school within each pair to receive treatment `Tr` and save RCT information. We also include a Bernoulli randomized treatment assignment `TrBern` (ignoring pairs).


set.seed(73)

rct <- poss_rct %>%
  select(CAMPUS) %>%
  mutate(match = fakematch) %>%
  filter(match %in% rct.ids) %>%
  mutate(match = factor(match,levels = rct.ids, labels = LETTERS[1:25])) %>%
  group_by(match) %>%
  mutate(k = row_number()) %>%
  pivot_wider(values_from = CAMPUS, names_from = k, names_prefix = "CAMPUS") %>%
  mutate(Tr1 = sample(0:1,1),
         Tr2 = 1-Tr1) %>%
  pivot_longer(!match,
               names_to = c(".value","k"),
               names_pattern = "(.*)(1|2)") %>%
  select(CAMPUS, match, k, Tr) %>%
  ungroup() %>%
  arrange(match, k) %>%
  left_join(select(covs_ms_noprep, CAMPUS, clust_size = CPETG08C_78))

rct$TrBern = sample(0:1, 50, replace = T)


## Rename Outcomes for Clarity
out_ms <- out_ms %>%
  rename(taks08 = outmA08) %>%
  select(CAMPUS, taks08, everything())


aux_dat = out_ms %>%
  filter(!(CAMPUS %in% rct$CAMPUS)) %>%
  left_join(covs_ms, by = "CAMPUS") %>%
  select(CAMPUS, GRDSPAN, taks08, starts_with("pre"),
         all_of(var.names.clean$var_clean), everything())


## SMALLER VIGNETTE VERSIONS

# Choose a subset of the variable names
num_vars <- 20 #number of variables to pull (may include missing indicators or not)

aux_dat_mod <- aux_dat %>%
  select(-ends_with("_34"), -ends_with("_34_mis")) %>%
  select(-ends_with("_45"), -ends_with("_45_mis")) %>%
  select(-ends_with("_56"), -ends_with("_56_mis")) %>%
  select(-ends_with("_78"), -ends_with("_78_mis")) %>%
  select(-starts_with("prem")) %>%
  select(-all_stud_n, -grade8_n) %>%
  select(-(starts_with("out") & ends_with("_na"))) %>%
  select(-starts_with("GRDSPAN")) %>%
  select(-exist34, -exist45, -exist56, -exist78) %>%
  select(-(starts_with("out") & ends_with("09")))



df_without_out <- aux_dat_mod %>%
  select(-starts_with("out"))

set.seed(29181)
rf <- randomForest(taks08 ~ . - CAMPUS, data = df_without_out)

var_scores <- caret::varImp(rf)

var_scores$var <- rownames(var_scores)

top_vars <- var_scores %>%
  arrange(desc(Overall)) %>%
  head(num_vars) %>%
  pull(var)

#get missing ones for those that have it and non-missing ones
small_var_names <- top_vars
for (i in 1:length(top_vars)) {
  if (grepl("_mis", top_vars[i], fixed = T)){
    non_mis_name <-str_sub(top_vars[i], start = 1, end = -5)
    if (non_mis_name %in% colnames(aux_dat)) {
      small_var_names <- c(small_var_names, non_mis_name)
    }
  }
  else {
    mis_name <- str_c(top_vars[i], "_mis")
    if (mis_name %in% colnames(aux_dat)) {
      small_var_names <- c(small_var_names, mis_name)
    }
  }
}

small_var_names <- unique(small_var_names)

#need to add background names and other outcome variables
outcome_names <- aux_dat %>%
  select(taks08, (starts_with("out") & !ends_with("na") & !ends_with("09"))) %>%
  colnames()


subset_var_names <- c(outcome_names, var.names.clean$var_clean, small_var_names)


# Update datasets

# version 1: no missing data, imputation already applied
schools_no_mis <- aux_dat %>%
  select(all_of(subset_var_names))

usethis::use_data(schools_no_mis, overwrite = TRUE)


# version 2: missing data, imputation reversed

for (col_name in colnames(schools_no_mis)) {
  if (str_detect(col_name, "_mis")) {
    missing_indices <- which(schools_no_mis[[col_name]] == 1)
    schools_no_mis[[str_remove(col_name, "_mis")]][missing_indices] <- NA
  }
}

# taks[C][A]_[g7][07]_[06]

# C, F, or M for campus, female or male
# A or M for all or math
# g7 or sum
#07 or 08 for the testing standard
#06 or 07 for the year of testing??

# TO DO - RENAME ANNOYING COLUMNS
schools <- select(schools_no_mis, -ends_with("_mis")) %>%
  select(-starts_with("outh")) %>% #remove high school outcomes
  rename(taksCM_sum08_07 = CA311PM07R_67,
         taksCM_sum07_07 = CA311TM07R_67,
         taksCA_sum07_07 = CA311TA07R_67,
         taksMM_sum08_07 = CM311PM07R_67,
         taksCM_g707_07 = CA007TM07R_67,
         taksCA_g707_07 = CA007TA07R_67,
         taksFA_sum07_07 = CF311TA07R_67,
         taksMM_sum07_07 = CM311TM07R_67,
         taksFM_sum08_07 = CF311PM07R_67,
         taksFM_sum07_07 = CF311TM07R_67,
         taksCA_sum08_07= CA311PA07R_67,
         taksMA_sum07_07 = CM311TA07R_67,
         taksMA_g707_07 = CM007TA07R_67,
         taksMM_g707_07= CM007TM07R_67,
         taksFA_sum08_07 = CF311PA07R_67,
         perc_campus_mobility = CPEMALLP_67,
         female_mobile_particip = CFMYA07R_67,
         taksMA_sum08_07 =CM311PA07R_67,
         taksCM_sum08_06 = CA311PM06R_67,
         taksCM_sum07_07_commended = CA311CM07R_67)

usethis::use_data(schools, overwrite = TRUE)
