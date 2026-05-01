############################################################
# LESION SITE SENSITIVITY ANALYSIS
#
# Tests whether fatigue-thalamic connectivity findings hold
# across lesion site subgroups in stroke participants.
#
# REQUIRES: Run Merge_and_demographics.R first so these
# objects exist in your environment:
#   - dat_model_gdss
#   - demo_subject
#   - out_dir
#
# PRODUCES:
#   - LESIONSITE_thalamic_full_stroke.csv
#   - LESIONSITE_thalamic_subcortical_only.csv
#   - LESIONSITE_thalamic_cortical_only.csv
#   - LESIONSITE_thalamic_site_moderation.csv
#   - LESIONSITE_summary_table.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(readr)
})

## ===================== CONFIG =====================
# Target circuits — must match feature_type/feature_name in dat_model_gdss
thalamic_circuits <- list(
  list(
    feature_type = "SEED_EDGE_FC",
    feature_name = "THalamus<->Insula",
    label        = "Thalamus-Insula",
    beta_main    = -0.0071,
    fdr_q        = 0.011
  ),
  list(
    feature_type = "SEED_TO_NET_FC",
    feature_name = "SEED:THalamus->VAN",
    label        = "Thalamus-VAN",
    beta_main    = -0.0043,
    fdr_q        = 0.012
  )
)

include_time <- TRUE   # match your main modeling script

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ===================== HELPERS =====================

ensure_pvalues <- function(td) {
  if (!("statistic" %in% names(td)) && ("t.value" %in% names(td)))
    td <- td %>% rename(statistic = t.value)
  if (!("p.value" %in% names(td)) && ("statistic" %in% names(td)))
    td <- td %>% mutate(p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE))
  td
}

p_fmt <- function(p) {
  sapply(p, function(x) {
    if (is.na(x))  return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  })
}

fit_thalamic_model <- function(df, circ, label_suffix) {
  dfp <- df %>%
    filter(
      feature_type == circ$feature_type,
      feature_name == circ$feature_name,
      !is.na(value),
      !is.na(behav_resid),
      !is.na(gdss_score),
      !is.na(age),
      !is.na(sex),
      !is.na(study_id)
    ) %>%
    { if (include_time) filter(., !is.na(time_mo)) else . } %>%
    droplevels()
  
  n_subj <- n_distinct(dfp$study_id)
  n_obs  <- nrow(dfp)
  
  if (n_obs < 8 || n_subj < 3) {
    cat(sprintf("  SKIPPED %s [%s]: only %d obs, %d subjects\n",
                circ$label, label_suffix, n_obs, n_subj))
    return(NULL)
  }
  
  form <- if (include_time) {
    as.formula("value ~ behav_resid * groupstroke + gdss_score + age + sex + time_mo + (1 | study_id)")
  } else {
    as.formula("value ~ behav_resid * groupstroke + gdss_score + age + sex + (1 | study_id)")
  }
  
  # For stroke-only subsets, drop group interaction
  if (n_distinct(as.character(dfp$group)) < 2) {
    form <- if (include_time) {
      as.formula("value ~ behav_resid + gdss_score + age + sex + time_mo + (1 | study_id)")
    } else {
      as.formula("value ~ behav_resid + gdss_score + age + sex + (1 | study_id)")
    }
  }
  
  if (nlevels(droplevels(dfp$sex)) < 2)
    form <- update(form, . ~ . - sex)
  
  m <- try(lmer(form, data = dfp, REML = FALSE, control = ctrl), silent = TRUE)
  if (inherits(m, "try-error")) {
    cat(sprintf("  FAILED  %s [%s]: %s\n",
                circ$label, label_suffix,
                conditionMessage(attr(m, "condition"))))
    return(NULL)
  }
  
  # Extract circuit metadata before mutate to avoid name collision
  circ_label      <- circ$label
  circ_beta_main  <- circ$beta_main
  circ_fdr_q      <- circ$fdr_q
  
  td <- broom.mixed::tidy(m, effects = "fixed") %>%
    ensure_pvalues() %>%
    mutate(
      circuit             = circ_label,
      subset              = label_suffix,
      n_obs               = n_obs,
      n_subjects          = n_subj,
      main_beta_published = circ_beta_main,
      main_fdr_published  = circ_fdr_q
    )
  
  cat(sprintf("  OK      %s [%s]: n=%d subjects=%d\n",
              circ_label, label_suffix, n_obs, n_subj))
  td
}

## ===================== 1) ATTACH LESION SITE =====================
cat("=== LESION SITE SENSITIVITY ANALYSIS ===\n\n")

# Attach lesion_site from demo_subject to dat_model_gdss
dat_site <- dat_model_gdss %>%
  left_join(
    demo_subject %>%
      select(study_id, lesion_site) %>%
      mutate(study_id = as.character(study_id)),
    by = "study_id"
  )

cat("Lesion site attached. Distribution in stroke rows:\n")
print(table(
  dat_site$lesion_site[as.character(dat_site$group) == "Stroke"],
  useNA = "ifany"
))
cat("\n")

## ===================== 2) RESIDUALIZE VITALITY =====================
# Residualize behav (vitality) for GDSS + covariates within each feature,
# consistent with the primary analysis approach

dat_site <- dat_site %>%
  group_by(feature_type, feature_name) %>%
  group_modify(~{
    d  <- .x
    ok <- is.finite(d$behav) & is.finite(d$gdss_score) &
      is.finite(d$age)   & !is.na(d$sex) & !is.na(d$group)
    if (include_time) ok <- ok & is.finite(d$time_mo)
    d$behav_resid <- NA_real_
    if (sum(ok) >= 10) {
      form_resid <- if (include_time) {
        behav ~ gdss_score + age + sex + group + time_mo
      } else {
        behav ~ gdss_score + age + sex + group
      }
      d$behav_resid[ok] <- resid(lm(form_resid, data = d[ok, , drop = FALSE]))
    }
    d
  }) %>%
  ungroup()

## ===================== 3) BUILD SUBSETS =====================

# Stroke only — full
dat_stroke_full <- dat_site %>%
  filter(as.character(group) == "Stroke")

# Stroke only — subcortical + cortico-subcortical
dat_subcort <- dat_site %>%
  filter(
    as.character(group) == "Stroke",
    lesion_site %in% c("Subcortical", "Cortico-subcortical",
                       "White matter only")
  )

# Stroke only — cortical only
dat_cortical <- dat_site %>%
  filter(
    as.character(group) == "Stroke",
    lesion_site == "Cortical"
  )

# Stroke only — posterior/other (brainstem, cerebellar, other)
dat_posterior <- dat_site %>%
  filter(
    as.character(group) == "Stroke",
    lesion_site %in% c("Brainstem", "Cerebellar", "Other")
  )

cat("Subset sizes (unique subjects):\n")
cat(sprintf("  Full stroke:            %d\n", n_distinct(dat_stroke_full$study_id)))
cat(sprintf("  Subcortical/WM:         %d\n", n_distinct(dat_subcort$study_id)))
cat(sprintf("  Cortical only:          %d\n", n_distinct(dat_cortical$study_id)))
cat(sprintf("  Posterior/other:        %d\n", n_distinct(dat_posterior$study_id)))
cat(sprintf("  Missing lesion site:    %d\n",
            n_distinct(dat_site$study_id[
              as.character(dat_site$group) == "Stroke" & is.na(dat_site$lesion_site)
            ])))
cat("\n")

## ===================== 4) RUN MODELS PER SUBSET =====================

subsets <- list(
  list(data = dat_stroke_full, label = "Full stroke (all sites)"),
  list(data = dat_subcort,     label = "Subcortical / cortico-subcortical / WM"),
  list(data = dat_cortical,    label = "Cortical only"),
  list(data = dat_posterior,   label = "Posterior / other")
)

all_results <- list()

for (ss in subsets) {
  cat(sprintf("--- Subset: %s ---\n", ss$label))
  for (circ in thalamic_circuits) {
    res <- fit_thalamic_model(ss$data, circ, ss$label)
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
  cat("\n")
}

results_df <- bind_rows(all_results)

## ===================== 5) LESION SITE MODERATION MODEL =====================
# Tests whether lesion_site_broad moderates the fatigue-connectivity relationship
# within stroke participants (Vitality_resid x lesion_site_broad interaction)

cat("--- Lesion site moderation models ---\n")

dat_stroke_mod <- dat_site %>%
  filter(as.character(group) == "Stroke", !is.na(lesion_site)) %>%
  mutate(
    lesion_site_broad = case_when(
      lesion_site %in% c("Subcortical", "White matter only") ~ "Subcortical",
      lesion_site == "Cortical"                              ~ "Cortical",
      lesion_site == "Cortico-subcortical"                   ~ "Mixed",
      lesion_site %in% c("Brainstem","Cerebellar","Other")   ~ "Posterior/Other",
      TRUE ~ NA_character_
    ),
    lesion_site_broad = factor(lesion_site_broad,
                               levels = c("Subcortical","Cortical",
                                          "Mixed","Posterior/Other"))
  ) %>%
  filter(!is.na(lesion_site_broad))

moderation_results <- list()

for (circ in thalamic_circuits) {
  dfp <- dat_stroke_mod %>%
    filter(
      feature_type == circ$feature_type,
      feature_name == circ$feature_name,
      !is.na(value), !is.na(behav_resid),
      !is.na(gdss_score), !is.na(age), !is.na(sex)
    ) %>%
    { if (include_time) filter(., !is.na(time_mo)) else . } %>%
    droplevels()
  
  n_obs  <- nrow(dfp)
  n_subj <- n_distinct(dfp$study_id)
  
  if (n_obs < 8 || n_subj < 5) {
    cat(sprintf("  SKIPPED moderation %s: n=%d\n", circ$label, n_obs))
    next
  }
  
  form_mod <- if (include_time) {
    as.formula("value ~ behav_resid * lesion_site_broad + gdss_score + age + sex + time_mo + (1 | study_id)")
  } else {
    as.formula("value ~ behav_resid * lesion_site_broad + gdss_score + age + sex + (1 | study_id)")
  }
  
  if (nlevels(droplevels(dfp$sex)) < 2)
    form_mod <- update(form_mod, . ~ . - sex)
  
  m_mod <- try(lmer(form_mod, data = dfp, REML = FALSE, control = ctrl),
               silent = TRUE)
  
  if (inherits(m_mod, "try-error")) {
    cat(sprintf("  FAILED moderation %s\n", circ$label))
    next
  }
  
  td_mod <- broom.mixed::tidy(m_mod, effects = "fixed") %>%
    ensure_pvalues() %>%
    mutate(
      circuit    = circ$label,
      subset     = "Lesion site moderation",
      n_obs      = n_obs,
      n_subjects = n_subj,
      main_beta_published = circ$beta_main,
      main_fdr_published  = circ$fdr_q
    )
  
  cat(sprintf("  OK moderation %s: n=%d subjects=%d\n",
              circ$label, n_obs, n_subj))
  moderation_results[[length(moderation_results) + 1]] <- td_mod
}

moderation_df <- bind_rows(moderation_results)
cat("\n")

## ===================== 6) SUMMARY TABLE =====================
# Extract the key behav_resid term from each subset for easy comparison

summary_table <- results_df %>%
  filter(term == "behav_resid") %>%
  select(circuit, subset, estimate, std.error, statistic, p.value,
         n_obs, n_subjects, main_beta_published, main_fdr_published) %>%
  mutate(
    p_label       = p_fmt(p.value),
    beta_change   = sprintf("%.4f vs %.4f (published)",
                            estimate, main_beta_published),
    direction_ok  = sign(estimate) == sign(main_beta_published)
  ) %>%
  arrange(circuit, subset)

cat("=== SENSITIVITY SUMMARY: behav_resid term across subsets ===\n")
print(summary_table %>%
        select(circuit, subset, estimate, p_label, n_subjects,
               direction_ok) %>%
        as.data.frame(),
      row.names = FALSE)
cat("\n")

# Moderation interaction terms
if (nrow(moderation_df) > 0) {
  cat("=== MODERATION: behav_resid x lesion_site_broad interactions ===\n")
  print(moderation_df %>%
          filter(str_detect(term, "behav_resid:")) %>%
          select(circuit, term, estimate, p.value, n_subjects) %>%
          mutate(p_label = p_fmt(p.value)) %>%
          as.data.frame(),
        row.names = FALSE)
  cat("\n")
}

## ===================== 7) SAVE OUTPUTS =====================

readr::write_csv(
  results_df,
  file.path(out_dir, "LESIONSITE_thalamic_all_subsets_full.csv")
)

readr::write_csv(
  summary_table,
  file.path(out_dir, "LESIONSITE_thalamic_summary_table.csv")
)

readr::write_csv(
  moderation_df,
  file.path(out_dir, "LESIONSITE_thalamic_site_moderation.csv")
)

# Separate files per subset for easy inspection
results_df %>%
  filter(subset == "Full stroke (all sites)") %>%
  readr::write_csv(file.path(out_dir, "LESIONSITE_thalamic_full_stroke.csv"))

results_df %>%
  filter(subset == "Subcortical / cortico-subcortical / WM") %>%
  readr::write_csv(file.path(out_dir, "LESIONSITE_thalamic_subcortical_only.csv"))

results_df %>%
  filter(subset == "Cortical only") %>%
  readr::write_csv(file.path(out_dir, "LESIONSITE_thalamic_cortical_only.csv"))

cat("\u2705 Lesion site sensitivity outputs saved:\n")
cat("   LESIONSITE_thalamic_all_subsets_full.csv\n")
cat("   LESIONSITE_thalamic_summary_table.csv\n")
cat("   LESIONSITE_thalamic_site_moderation.csv\n")
cat("   LESIONSITE_thalamic_full_stroke.csv\n")
cat("   LESIONSITE_thalamic_subcortical_only.csv\n")
cat("   LESIONSITE_thalamic_cortical_only.csv\n")

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("Key questions to answer from summary_table:\n")
cat("  1. direction_ok == TRUE in all subsets?\n")
cat("     -> Effect direction consistent with published finding\n")
cat("  2. Significant (p < 0.05) in subcortical/cortico-subcortical subset?\n")
cat("     -> Effect holds where direct thalamic pathway disruption is most plausible\n")
cat("  3. Significant in cortical-only subset?\n")
cat("     -> Effect holds even without direct thalamic lesion (strongest network argument)\n")
cat("  4. Moderation interaction significant?\n")
cat("     -> Lesion location moderates the fatigue-connectivity relationship\n")
cat("     -> If NOT significant: effect is lesion-site independent (strengthens claim)\n")