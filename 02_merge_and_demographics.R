#=============================================================================
  #  Script:  merge_and_demographics.R
  #
  #  Purpose: Merges behavioural/REDCap data with imaging feature data,
  #           extracts demographic characteristics, computes group comparison
  #           statistics, and exports analysis-ready long-format CSVs.
  #
  #  Paper:   "Depression-Independent Thalamocortical Correlates of
  #            Post-Stroke Fatigue"
  #
  #  Outputs:
  #    - merged_imaging_behaviour_long_ALL.csv       (all subjects)
  #    - merged_imaging_behaviour_long_MODEL.csv     (vitality complete)
  #    - merged_imaging_behaviour_long_MODEL_BEHAV_GDSS.csv (vitality+GDSS)
  #    - DEMOGRAPHICS_table_manuscript.csv
  #    - DEMOGRAPHICS_subject_level.csv
  #    - DEMOGRAPHICS_timepoint_attendance.csv
  #    - DEMOGRAPHICS_vitality_depression_correlation.csv
  #
  #  Requires:
  #    readxl, readr, dplyr, tidyr, stringr, purrr, broom, lme4, broom.mixed
  #
  #  Input files:
  #    - behav_path : REDCap export (.xlsx) containing demographics,
  #                   SF-36 Vitality (sf36_en), and GDSS scores
  #    - img_path   : Long-format imaging features CSV with columns:
  #                   study_id, event_name, feature_type, feature_name, value
  #
  #  Usage:
  #    1. Set the three paths in the CONFIG section below
  #    2. Verify the RECODE BLOCKS match your REDCap numeric codes
  #    3. Source the entire script
  #
  #  Note: Participant data are not publicly available due to confidentiality.
  #        Code is shared for methodological transparency.
  #        Contact the corresponding author to request data access.
  # =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(broom)
  library(lme4)
  library(broom.mixed)
})

## ===================== CONFIG =====================
behav_path <- "/Users/.../FCS_Demographics_and_Behavior_Dec_19_2017_arm1.xlsx"
img_path   <- "/Users/.../regionFC_and_ChaCo_alff_falff_features_long.csv"
out_dir    <- "/Users/..."

behav_var <- "sf36_en"
gdss_var  <- "gdss_score"
id_col    <- "study_id"
event_col <- "redcap_event_name"
group_col <- "subj_type"
sex_col   <- "gender"
age_col   <- "age"

feature_types_keep <- c(
  "WB_FC","NET_BETWEEN_FC","NET_WITHIN_FC","FC_NODE_STRENGTH",
  "DLPFC_TO_NET_FC","DLPFC_TO_ROI_FC",
  "SEED_EDGE_FC","SEED_TO_NET_FC",
  "ALFF_ROI","FALFF_ROI","ALFF_NET_MEAN","FALFF_NET_MEAN",
  "CHACO_NET_MEAN","CHACO_NET_WITHIN_MEAN","CHACO_NET_BETWEEN_MEAN"
)

## ===================== HELPERS (unchanged) =====================

standardize_id <- function(x) {
  s <- x %>%
    as.character() %>%
    str_trim() %>%
    tolower() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("-", "_")
  s <- s %>% str_replace("^fcs_?0*([0-9]+)$", "fcs_\\1")
  s <- ifelse(
    str_detect(s, "^fcs_\\d+$"),
    paste0("fcs_", str_pad(str_extract(s, "(?<=fcs_)\\d+"), 3, pad = "0")),
    s
  )
  s <- ifelse(
    str_detect(s, "^amc_\\d+_fcs$"),
    paste0("amc_", str_pad(str_extract(s, "(?<=amc_)\\d+"), 3, pad = "0"), "_fcs"),
    s
  )
  s
}

clean_missing <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x[x %in% c("", "NA", "N/A", "NaN", "nan")] <- NA
  x
}

pick_first_nonmissing <- function(x) {
  x <- clean_missing(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

harmonize_event <- function(ev) {
  e <- ev %>% as.character() %>% str_trim() %>% tolower()
  dplyr::case_when(
    e == "acute_arm_1"                ~ "acute_arm_1",
    e == "3month_arm_1"               ~ "3month_arm_1",
    e == "1year_arm_1"                ~ "1year_arm_1",
    str_detect(e, "^visit_1_arm_2")   ~ "visit_1_arm_2",
    str_detect(e, "^visit_2_arm_2")   ~ "visit_2_arm_2",
    TRUE ~ NA_character_
  )
}

event_to_time <- function(ev) {
  case_when(
    ev == "acute_arm_1"               ~ 0,
    ev == "3month_arm_1"              ~ 3,
    ev == "1year_arm_1"               ~ 12,
    str_detect(ev, "^visit_1")        ~ 1,
    str_detect(ev, "^visit_2")        ~ 2,
    TRUE ~ NA_real_
  )
}

## ===================== DEMOGRAPHICS HELPERS =====================

fmt_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
fmt_median_iqr <- function(x) {
  sprintf("%.1f [%.1f\u2013%.1f]",
          median(x,             na.rm = TRUE),
          quantile(x, 0.25,     na.rm = TRUE),
          quantile(x, 0.75,     na.rm = TRUE))
}
fmt_n_pct <- function(condition, total) {
  n <- sum(condition, na.rm = TRUE)
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}
p_fmt <- function(p) {
  sapply(p, function(x) {
    if (is.na(x))  return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  })
}

## ===================== 1) LOAD BEHAVIOR =====================
behav_raw <- readxl::read_excel(behav_path)

stopifnot(all(c(id_col, event_col, group_col, sex_col, age_col,
                behav_var, gdss_var) %in% names(behav_raw)))

# Print unique values of key coding columns so you can verify recodes below
cat("=== KEY COLUMN VALUES (verify recodes below match) ===\n")
for (col in c(group_col, sex_col, "race", "ethnicity", "handed",
              "education", "nihss_hospital", "lesion_side", "lesion_site",
              "lesion_type", "tpa", "thrombectomy")) {
  if (col %in% names(behav_raw)) {
    vals <- sort(unique(behav_raw[[col]]))
    cat(sprintf("  %-20s: %s\n", col,
                paste(head(vals, 15), collapse = ", ")))
  }
}
cat("\n")

behav <- behav_raw %>%
  mutate(
    study_id       = standardize_id(.data[[id_col]]),
    event_name_raw = as.character(.data[[event_col]]),
    event_name     = harmonize_event(event_name_raw),
    group_raw      = .data[[group_col]],
    sex_raw        = .data[[sex_col]],
    age_raw        = .data[[age_col]],
    behav_raw      = clean_missing(.data[[behav_var]]),
    gdss_raw       = clean_missing(.data[[gdss_var]])
  )

cat("Behavior rows:", nrow(behav), "\n")
cat("Behavior harmonized events:\n")
print(sort(unique(behav$event_name)))
cat("Behavior subjects:", n_distinct(behav$study_id), "\n\n")

## ===================== 2) DEMOGRAPHICS EXTRACTION =====================
# Pull from basic_subject_info row first; fall back to any row per subject
# RECODE BLOCKS: edit numeric codes below to match the printout above

demo_subject <- behav %>%
  mutate(is_basic = str_detect(
    tolower(str_trim(event_name_raw)), "^basic_subject_info")) %>%
  group_by(study_id) %>%
  summarise(
    # ---- group ----
    group_raw_val = {
      g <- pick_first_nonmissing(group_raw[is_basic])
      if (is.na(g)) pick_first_nonmissing(group_raw) else g
    },
    # ---- sex ----
    sex_raw_val = {
      s <- pick_first_nonmissing(sex_raw[is_basic])
      if (is.na(s)) pick_first_nonmissing(sex_raw) else s
    },
    # ---- age ----
    age = suppressWarnings(as.numeric({
      a <- pick_first_nonmissing(age_raw[is_basic])
      if (is.na(a)) pick_first_nonmissing(age_raw) else a
    })),
    # ---- additional demographics (basic_subject_info row only) ----
    race_raw       = pick_first_nonmissing(if ("race"       %in% names(cur_data())) .data$race[is_basic]       else rep(NA, sum(is_basic))),
    ethnicity_raw  = pick_first_nonmissing(if ("ethnicity"  %in% names(cur_data())) .data$ethnicity[is_basic]  else rep(NA, sum(is_basic))),
    handed_raw     = pick_first_nonmissing(if ("handed"     %in% names(cur_data())) .data$handed[is_basic]     else rep(NA, sum(is_basic))),
    education_raw  = pick_first_nonmissing(if ("education"  %in% names(cur_data())) .data$education[is_basic]  else rep(NA, sum(is_basic))),
    # ---- stroke clinical (acute row preferred) ----
    nihss          = suppressWarnings(as.numeric(
      pick_first_nonmissing(if ("nihss_hospital" %in% names(cur_data())) .data$nihss_hospital else rep(NA, n()))
    )),
    lesion_size    = suppressWarnings(as.numeric(
      pick_first_nonmissing(if ("lesion_size"    %in% names(cur_data())) .data$lesion_size    else rep(NA, n()))
    )),
    lesion_side_raw = pick_first_nonmissing(
      if ("lesion_side" %in% names(cur_data())) .data$lesion_side else rep(NA, n())
    ),
    lesion_type_raw = pick_first_nonmissing(
      if ("lesion_type" %in% names(cur_data())) .data$lesion_type else rep(NA, n())
    ),
    tpa_raw        = pick_first_nonmissing(
      if ("tpa"         %in% names(cur_data())) .data$tpa         else rep(NA, n())
    ),
    thrombectomy_raw = pick_first_nonmissing(
      if ("thrombectomy" %in% names(cur_data())) .data$thrombectomy else rep(NA, n())
    ),
    lesion_site_raw = pick_first_nonmissing(
      if ("lesion_site" %in% names(cur_data())) .data$lesion_site else rep(NA, n())
    ),
    # ---- store event arm for group fallback ----
    event_arm = pick_first_nonmissing(event_name_raw),
    .groups = "drop"
  ) %>%
  mutate(
    # ================================================================
    # RECODE BLOCKS — edit these to match the printout above
    # ================================================================
    
    # GROUP: subj_type 0 = stroke, 1 = control
    # Fallback: arm_2 events = control, arm_1 events = stroke
    # (controls are in REDCap arm_2 and may have NA for subj_type)
    group = case_when(
      as.character(group_raw_val) == "0"                    ~ "Stroke",
      as.character(group_raw_val) == "1"                    ~ "Control",
      TRUE ~ NA_character_
    ),
    
    # SEX: gender 0 = female, 1 = male (adjust if different)
    sex = case_when(
      as.character(sex_raw_val) == "1" ~ "Female",
      as.character(sex_raw_val) == "0" ~ "Male",
      TRUE ~ NA_character_
    ),
    
    # RACE (adjust numeric codes to match your REDCap)
    race = case_when(
      as.character(race_raw) == "4" ~ "White",
      as.character(race_raw) == "3" ~ "Black or African American",
      as.character(race_raw) == "1" ~ "Asian",
      as.character(race_raw) == "0" ~ "American Indian or Alaska Native",
      as.character(race_raw) == "2" ~ "Native Hawaiian or Pacific Islander",
      TRUE ~ NA_character_
    ),
    
    # ETHNICITY: 0 = not Hispanic, 1 = Hispanic (adjust if different)
    ethnicity = case_when(
      as.character(ethnicity_raw) == "1" ~ "Not Hispanic or Latino",
      as.character(ethnicity_raw) == "0" ~ "Hispanic or Latino",
      TRUE ~ NA_character_
    ),
    
    # HANDEDNESS: 0 = right, 1 = left, 2 = ambidextrous (adjust if different)
    handedness = case_when(
      as.character(handed_raw) == "0" ~ "Left",
      as.character(handed_raw) == "1" ~ "Right",
      as.character(handed_raw) == "2" ~ "Ambidextrous",
      TRUE ~ NA_character_
    ),
    
    education = suppressWarnings(as.numeric(education_raw)),
    
    # LESION SIDE
    lesion_side = case_when(
      as.character(lesion_side_raw) == "0" ~ "Left",
      as.character(lesion_side_raw) == "1" ~ "Right",
      TRUE ~ NA_character_
    ),
    
    # LESION TYPE
    lesion_type = case_when(
      as.character(lesion_type_raw) == "0" ~ "Ischemic",
      as.character(lesion_type_raw) == "1" ~ "Hemorrhagic",
      TRUE ~ NA_character_
    ),
    
    # LESION SITE
    lesion_site = case_when(
      as.character(lesion_site_raw) == "0" ~ "Subcortical",
      as.character(lesion_site_raw) == "1" ~ "Cortical",
      as.character(lesion_site_raw) == "2" ~ "Cortico-subcortical",
      as.character(lesion_site_raw) == "3" ~ "White matter only",
      as.character(lesion_site_raw) == "4" ~ "Brainstem",
      as.character(lesion_site_raw) == "5" ~ "Cerebellar",
      as.character(lesion_site_raw) == "6" ~ "Other",
      TRUE ~ NA_character_
    ),
    
    tpa          = as.numeric(as.character(tpa_raw)) == 1,
    thrombectomy = as.numeric(as.character(thrombectomy_raw)) == 1
  )

cat("demo_subject N =", nrow(demo_subject), "\n")
cat("Group counts (including arm fallback):\n")
print(table(demo_subject$group, useNA = "ifany"))
cat("Sex counts:\n")
print(table(demo_subject$sex, useNA = "ifany"))
cat("Lesion site counts (stroke only):\n")
print(table(demo_subject$lesion_site[demo_subject$group == "Stroke"], useNA = "ifany"))
cat("\n")

## ===================== 3) PER-EVENT BEHAVIOR =====================
behav_event <- behav %>%
  filter(!is.na(event_name)) %>%
  group_by(study_id, event_name) %>%
  summarise(
    behav      = pick_first_nonmissing(behav_raw),
    gdss_score = pick_first_nonmissing(gdss_raw),
    .groups    = "drop"
  )

behav_clean <- behav_event %>%
  left_join(demo_subject %>%
              select(study_id, group, sex, age),
            by = "study_id") %>%
  mutate(
    time_mo    = event_to_time(event_name),
    group      = factor(group),
    sex        = factor(sex),
    age        = suppressWarnings(as.numeric(age)),
    behav      = suppressWarnings(as.numeric(behav)),
    gdss_score = suppressWarnings(as.numeric(gdss_score))
  )

cat("After cleaning behavior:\n")
cat("  rows:", nrow(behav_clean), "\n")
cat("  missing behav:", sum(is.na(behav_clean$behav)), "\n")
cat("  missing gdss :", sum(is.na(behav_clean$gdss_score)), "\n")

## ===================== 4) SUBJECT-LEVEL BEHAVIOURAL SUMMARY =====================
subj_beh <- behav_clean %>%
  filter(!is.na(behav) | !is.na(gdss_score)) %>%
  group_by(study_id) %>%
  summarise(
    vitality_mean = ifelse(all(is.na(behav)),      NA_real_, mean(behav,      na.rm = TRUE)),
    vitality_base = suppressWarnings(behav[which.min(replace(time_mo, is.na(time_mo), Inf))][1]),
    gdss_mean     = ifelse(all(is.na(gdss_score)), NA_real_, mean(gdss_score, na.rm = TRUE)),
    gdss_base     = suppressWarnings(gdss_score[which.min(replace(time_mo, is.na(time_mo), Inf))][1]),
    n_timepoints  = n_distinct(time_mo[!is.na(behav) & !is.na(gdss_score)]),
    timepoints_present = paste(sort(unique(na.omit(time_mo[!is.na(behav)]))),
                               collapse = ", "),
    .groups = "drop"
  )

# Full subject demographics + behavioural summary
subj_all <- demo_subject %>%
  left_join(subj_beh, by = "study_id")

cat("\nSubject-level dataset N =", nrow(subj_all), "\n")
cat("With vitality_mean:    ", sum(!is.na(subj_all$vitality_mean)), "\n")
cat("With gdss_mean:        ", sum(!is.na(subj_all$gdss_mean)), "\n\n")

## ===================== 5) TIMEPOINT ATTENDANCE =====================
tp_attend <- behav_clean %>%
  filter(!is.na(behav), !is.na(gdss_score), !is.na(group)) %>%
  distinct(study_id, group, time_mo) %>%
  group_by(group, time_mo) %>%
  summarise(n = n_distinct(study_id), .groups = "drop") %>%
  arrange(group, time_mo)

cat("=== TIMEPOINT ATTENDANCE (subjects with vitality + GDSS) ===\n")
print(tp_attend)
cat("\n")

## ===================== 6) VITALITY-DEPRESSION CORRELATIONS =====================
cor_data <- behav_clean %>%
  filter(!is.na(behav), !is.na(gdss_score),
         !is.na(group), group %in% c("Stroke", "Control")) %>%
  distinct(study_id, time_mo, behav, gdss_score, group)

cat(sprintf("Correlation dataset: %d observations, %d subjects\n",
            nrow(cor_data), n_distinct(cor_data$study_id)))

# Observation-level
cor_obs_p <- cor.test(cor_data$behav, cor_data$gdss_score, method = "pearson")
cor_obs_s <- suppressWarnings(
  cor.test(cor_data$behav, cor_data$gdss_score, method = "spearman"))

# By group
cor_by_group <- cor_data %>%
  group_by(group) %>%
  summarise(
    n          = n(),
    pearson_r  = cor(behav, gdss_score, use = "complete.obs", method = "pearson"),
    pearson_p  = cor.test(behav, gdss_score, method = "pearson")$p.value,
    spearman_r = cor(behav, gdss_score, use = "complete.obs", method = "spearman"),
    spearman_p = suppressWarnings(
      cor.test(behav, gdss_score, method = "spearman")$p.value),
    .groups = "drop"
  )

# Subject-level
cor_subj <- subj_all %>%
  filter(!is.na(vitality_mean), !is.na(gdss_mean),
         !is.na(group), as.character(group) %in% c("Stroke", "Control"))

cat(sprintf("Subject-level correlation n = %d\n", nrow(cor_subj)))
cor_subj_p <- cor.test(cor_subj$vitality_mean, cor_subj$gdss_mean, method = "pearson")
cor_subj_s <- suppressWarnings(
  cor.test(cor_subj$vitality_mean, cor_subj$gdss_mean, method = "spearman"))

cat("=== VITALITY-DEPRESSION CORRELATIONS ===\n")
cat("Observation-level:\n")
cat(sprintf("  Pearson  r = %.3f, p = %s\n",
            cor_obs_p$estimate, p_fmt(cor_obs_p$p.value)))
cat(sprintf("  Spearman r = %.3f, p = %s\n",
            cor_obs_s$estimate, p_fmt(cor_obs_s$p.value)))
cat("By group:\n")
print(cor_by_group %>%
        mutate(pearson_p  = p_fmt(pearson_p),
               spearman_p = p_fmt(spearman_p)))
cat(sprintf("Subject-level (n = %d):\n", nrow(cor_subj)))
cat(sprintf("  Pearson  r = %.3f, p = %s\n",
            cor_subj_p$estimate, p_fmt(cor_subj_p$p.value)))
cat(sprintf("  Spearman r = %.3f, p = %s\n\n",
            cor_subj_s$estimate, p_fmt(cor_subj_s$p.value)))

## ===================== 7) GROUP COMPARISON TESTS =====================
grp_sub <- subj_all %>%
  filter(!is.na(group), as.character(group) %in% c("Stroke", "Control")) %>%
  mutate(group = factor(as.character(group), levels = c("Stroke", "Control"))) %>%
  droplevels()

cat("Group levels in grp_sub:", levels(grp_sub$group), "\n")
cat("Group counts:\n")
print(table(grp_sub$group))
cat("\n")

grp_age <- wilcox.test(age           ~ group, data = grp_sub)
grp_vit <- wilcox.test(vitality_mean ~ group, data = grp_sub)
grp_gds <- wilcox.test(gdss_mean     ~ group, data = grp_sub)
grp_sex <- chisq.test(table(
  grp_sub$sex[!is.na(grp_sub$sex)],
  grp_sub$group[!is.na(grp_sub$sex)]
))

cat("=== GROUP COMPARISON TESTS ===\n")
cat(sprintf("Age:      Wilcoxon p = %s\n", p_fmt(grp_age$p.value)))
cat(sprintf("Sex:      Chi-sq   p = %s\n", p_fmt(grp_sex$p.value)))
cat(sprintf("Vitality: Wilcoxon p = %s\n", p_fmt(grp_vit$p.value)))
cat(sprintf("GDSS:     Wilcoxon p = %s\n\n", p_fmt(grp_gds$p.value)))

## ===================== 8) DEMOGRAPHICS TABLE =====================

make_demo_block <- function(df, label) {
  N <- nrow(df)
  tibble(
    Variable = c(
      "N",
      "Age mean (SD)", "Age median [IQR]",
      "Female n (%)", "Male n (%)",
      "White n (%)", "Black or African American n (%)",
      "Hispanic or Latino n (%)",
      "Right-handed n (%)",
      "Education years mean (SD)",
      "Vitality mean (SD)", "Vitality median [IQR]",
      "GDSS mean (SD)", "GDSS median [IQR]",
      "1 timepoint n (%)", "2 timepoints n (%)", "3 timepoints n (%)"
    ),
    !!label := c(
      as.character(N),
      fmt_mean_sd(df$age),               fmt_median_iqr(df$age),
      fmt_n_pct(df$sex == "Female", N),  fmt_n_pct(df$sex == "Male", N),
      fmt_n_pct(df$race == "White", N),
      fmt_n_pct(df$race == "Black or African American", N),
      fmt_n_pct(df$ethnicity == "Hispanic or Latino", N),
      fmt_n_pct(df$handedness == "Right", N),
      fmt_mean_sd(df$education),
      fmt_mean_sd(df$vitality_mean),     fmt_median_iqr(df$vitality_mean),
      fmt_mean_sd(df$gdss_mean),         fmt_median_iqr(df$gdss_mean),
      fmt_n_pct(df$n_timepoints == 1, N),
      fmt_n_pct(df$n_timepoints == 2, N),
      fmt_n_pct(df$n_timepoints == 3, N)
    )
  )
}

make_stroke_extras <- function(df, label) {
  N <- nrow(df)
  tibble(
    Variable = c(
      "NIHSS median [IQR]",
      "Lesion size cm3 median [IQR]",
      "Left hemisphere n (%)", "Right hemisphere n (%)", "Bilateral n (%)",
      "Ischemic n (%)", "Hemorrhagic n (%)",
      "tPA n (%)", "Thrombectomy n (%)",
      "Lesion site: Subcortical n (%)",
      "Lesion site: Cortical n (%)",
      "Lesion site: Cortico-subcortical n (%)",
      "Lesion site: White matter only n (%)",
      "Lesion site: Brainstem n (%)",
      "Lesion site: Cerebellar n (%)",
      "Lesion site: Other n (%)"
    ),
    !!label := c(
      fmt_median_iqr(df$nihss),
      fmt_median_iqr(df$lesion_size),
      fmt_n_pct(df$lesion_side == "Left",            N),
      fmt_n_pct(df$lesion_side == "Right",           N),
      fmt_n_pct(df$lesion_side == "Bilateral",       N),
      fmt_n_pct(df$lesion_type == "Ischemic",        N),
      fmt_n_pct(df$lesion_type == "Hemorrhagic",     N),
      fmt_n_pct(df$tpa          == TRUE,             N),
      fmt_n_pct(df$thrombectomy == TRUE,             N),
      fmt_n_pct(df$lesion_site == "Subcortical",          N),
      fmt_n_pct(df$lesion_site == "Cortical",             N),
      fmt_n_pct(df$lesion_site == "Cortico-subcortical",  N),
      fmt_n_pct(df$lesion_site == "White matter only",    N),
      fmt_n_pct(df$lesion_site == "Brainstem",            N),
      fmt_n_pct(df$lesion_site == "Cerebellar",           N),
      fmt_n_pct(df$lesion_site == "Other",                N)
    )
  )
}

stroke_df  <- subj_all %>% filter(!is.na(group), as.character(group) == "Stroke")
control_df <- subj_all %>% filter(!is.na(group), as.character(group) == "Control")
total_df   <- subj_all %>% filter(!is.na(group), as.character(group) %in% c("Stroke","Control"))

tbl_stroke       <- make_demo_block(stroke_df,  "Stroke")
tbl_control      <- make_demo_block(control_df, "Control")
tbl_total        <- make_demo_block(total_df,   "Total")
tbl_stroke_extra <- make_stroke_extras(stroke_df, "Stroke") %>%
  mutate(Control = "\u2014", Total = "\u2014")

demo_final <- tbl_stroke %>%
  left_join(tbl_control, by = "Variable") %>%
  left_join(tbl_total,   by = "Variable") %>%
  bind_rows(tbl_stroke_extra)

cat("=== DEMOGRAPHICS TABLE ===\n")
print(demo_final, n = Inf)

## ===================== 9) SAVE DEMOGRAPHICS OUTPUTS =====================
readr::write_csv(
  subj_all %>% select(-ends_with("_raw"), -ends_with("_raw_val")),
  file.path(out_dir, "DEMOGRAPHICS_subject_level.csv")
)

readr::write_csv(
  demo_final,
  file.path(out_dir, "DEMOGRAPHICS_table_manuscript.csv")
)

readr::write_csv(
  tp_attend,
  file.path(out_dir, "DEMOGRAPHICS_timepoint_attendance.csv")
)

readr::write_csv(
  tibble(
    level   = c("observation_pooled","observation_pooled",
                "subject_mean",      "subject_mean"),
    method  = c("pearson","spearman","pearson","spearman"),
    r       = c(cor_obs_p$estimate, cor_obs_s$estimate,
                cor_subj_p$estimate, cor_subj_s$estimate),
    p       = c(cor_obs_p$p.value, cor_obs_s$p.value,
                cor_subj_p$p.value, cor_subj_s$p.value),
    p_label = p_fmt(c(cor_obs_p$p.value, cor_obs_s$p.value,
                      cor_subj_p$p.value, cor_subj_s$p.value))
  ),
  file.path(out_dir, "DEMOGRAPHICS_vitality_depression_correlation.csv")
)

cat("\n\u2705 Demographics outputs saved.\n\n")

## ===================== 10) IMAGING PIPELINE (unchanged from original) =====================

norm_ft <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- trimws(x)
  toupper(x)
}
norm_fn <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  trimws(x)
}

keep_ft <- norm_ft(c(
  "WB_FC","NET_BETWEEN_FC","NET_WITHIN_FC","FC_NODE_STRENGTH",
  "DLPFC_TO_NET_FC","DLPFC_TO_ROI_FC",
  "SEED_EDGE_FC","SEED_TO_NET_FC",
  "ALFF_ROI","FALFF_ROI","ALFF_NET_MEAN","FALFF_NET_MEAN",
  "CHACO_NET_MEAN","CHACO_NET_WITHIN_MEAN","CHACO_NET_BETWEEN_MEAN"
))

rm(list = intersect("img", ls()))
img <- readr::read_csv(img_path, show_col_types = FALSE)

stopifnot(all(c("study_id","event_name","feature_type",
                "feature_name","value") %in% names(img)))

img <- img %>%
  mutate(
    study_id     = standardize_id(study_id),
    event_name   = as.character(event_name),
    feature_type = norm_ft(feature_type),
    feature_name = norm_fn(feature_name),
    value        = suppressWarnings(as.numeric(value))
  )

cat("\n=== IMAGING DIAGNOSTICS ===\n")
cat("Raw imaging rows:", nrow(img), "\n")
present_ft    <- sort(unique(img$feature_type))
missing_ft    <- setdiff(keep_ft, present_ft)
cat("Feature types present:\n"); print(present_ft)
cat("Expected but missing:\n");  print(missing_ft)

img <- img %>% filter(feature_type %in% keep_ft)
cat("Rows after keep filter:", nrow(img), "\n\n")

dup_img <- img %>%
  count(study_id, event_name, feature_type, feature_name, name = "n") %>%
  filter(n > 1)
cat("Duplicate imaging keys:", nrow(dup_img), "\n")

dat <- img %>%
  left_join(behav_clean, by = c("study_id", "event_name"))

cat("\nMerged rows:", nrow(dat), "\n")
cat("Missing behav:", sum(is.na(dat$behav)),      "\n")
cat("Missing gdss :", sum(is.na(dat$gdss_score)), "\n")
cat("Missing age  :", sum(is.na(dat$age)),        "\n")
cat("Missing group:", sum(is.na(dat$group)),      "\n\n")

dat0 <- dat %>%
  filter(!is.na(value)) %>%
  mutate(
    subject    = factor(study_id),
    event_name = factor(event_name),
    seed_name  = case_when(
      feature_type %in% c("SEED_TO_NET_FC","SEED_TO_ROI_FC") ~
        str_match(feature_name, "^SEED:([^->]+)")[,2],
      feature_type == "SEED_EDGE_FC" ~
        str_match(feature_name, "^([^<]+)<->([^<]+)$")[,2],
      TRUE ~ NA_character_
    ),
    model_family = case_when(
      feature_type %in% c("WB_FC","NET_BETWEEN_FC","NET_WITHIN_FC") ~ "FC_NETWORK",
      str_detect(feature_type, "^DLPFC_")                           ~ "FC_DLPFC",
      feature_type %in% c("SEED_TO_NET_FC","SEED_TO_ROI_FC",
                          "SEED_EDGE_FC")                          ~ "FC_SEED",
      feature_type == "FC_NODE_STRENGTH"                            ~ "FC_NODE_STRENGTH",
      feature_type == "ALFF_NET_MEAN"                               ~ "ALFF_NETWORK",
      feature_type == "FALFF_NET_MEAN"                              ~ "FALFF_NETWORK",
      feature_type == "ALFF_ROI"                                    ~ "ALFF_ROI",
      feature_type == "FALFF_ROI"                                   ~ "FALFF_ROI",
      str_detect(feature_type, "^CHACO")                            ~ "CHACO",
      TRUE ~ "OTHER"
    )
  )

dat_all        <- dat0 %>% filter(!is.na(age), !is.na(sex), !is.na(group))
dat_model      <- dat_all %>% filter(!is.na(behav))
dat_model_gdss <- dat_all %>% filter(!is.na(behav), !is.na(gdss_score))

cat("Rows dat_all:        ", nrow(dat_all),        "\n")
cat("Rows dat_model:      ", nrow(dat_model),      "\n")
cat("Rows dat_model_gdss: ", nrow(dat_model_gdss), "\n\n")

## ===================== 11) EXPORT MERGED CSVs =====================

write_merged <- function(df, fname) {
  df %>%
    select(study_id, event_name, time_mo,
           feature_type, feature_name, model_family, seed_name,
           value, behav, gdss_score, age, sex, group) %>%
    arrange(study_id, event_name, feature_type, feature_name) %>%
    readr::write_csv(file.path(out_dir, fname))
}

write_merged(dat_all,        "merged_imaging_behaviour_long_ALL.csv")
write_merged(dat_model,      "merged_imaging_behaviour_long_MODEL.csv")
write_merged(dat_model_gdss, "merged_imaging_behaviour_long_MODEL_BEHAV_GDSS.csv")

cat("\u2705 Merged CSVs saved.\n")
cat("   merged_imaging_behaviour_long_ALL.csv\n")
cat("   merged_imaging_behaviour_long_MODEL.csv\n")
cat("   merged_imaging_behaviour_long_MODEL_BEHAV_GDSS.csv\n")