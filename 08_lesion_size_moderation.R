############################################################
# LESION SIZE MODERATION ANALYSIS
#
# Tests whether lesion size moderates the fatigue-thalamic
# connectivity relationship in stroke participants.
#
# Analogous to ChaCo moderation — continuous moderator,
# stroke-only, residualized vitality as behavioral predictor.
#
# REQUIRES: Run Merge_and_demographics.R first so these
# objects exist in your environment:
#   - dat_model_gdss
#   - demo_subject
#   - out_dir
#
# PRODUCES:
#   - LESIONSIZE_thalamic_moderation_full.csv
#   - LESIONSIZE_thalamic_moderation_summary.csv
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

## ===================== CONFIG =====================
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

include_time <- TRUE
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

## ===================== 1) ATTACH LESION SIZE + RESIDUALIZE =====================
cat("=== LESION SIZE MODERATION ANALYSIS ===\n\n")

# Attach lesion_size from demo_subject
dat_ls <- dat_model_gdss %>%
  left_join(
    demo_subject %>%
      select(study_id, lesion_size) %>%
      mutate(study_id = as.character(study_id)),
    by = "study_id"
  )

# Stroke only, with lesion size present
dat_stroke_ls <- dat_ls %>%
  filter(as.character(group) == "Stroke",
         !is.na(lesion_size),
         lesion_size > 0)

cat(sprintf("Stroke subjects with lesion size: %d\n",
            n_distinct(dat_stroke_ls$study_id)))
cat(sprintf("Lesion size range: %.1f to %.1f cm3\n",
            min(dat_stroke_ls$lesion_size, na.rm = TRUE),
            max(dat_stroke_ls$lesion_size, na.rm = TRUE)))
cat(sprintf("Lesion size median [IQR]: %.1f [%.1f, %.1f] cm3\n\n",
            median(dat_stroke_ls$lesion_size, na.rm = TRUE),
            quantile(dat_stroke_ls$lesion_size, 0.25, na.rm = TRUE),
            quantile(dat_stroke_ls$lesion_size, 0.75, na.rm = TRUE)))

# Note: lesion size is right-skewed — log-transform for moderation
# This makes the interaction term more interpretable and reduces
# leverage of very large lesions
dat_stroke_ls <- dat_stroke_ls %>%
  mutate(
    lesion_size_log = log(lesion_size),
    lesion_size_z   = as.numeric(scale(lesion_size_log))
  )

cat("Lesion size log-transformed and z-scored for moderation models.\n\n")

# Residualize vitality within each feature
dat_stroke_ls <- dat_stroke_ls %>%
  group_by(feature_type, feature_name) %>%
  group_modify(~{
    d  <- .x
    ok <- is.finite(d$behav) & is.finite(d$gdss_score) &
      is.finite(d$age)   & !is.na(d$sex)
    if (include_time) ok <- ok & is.finite(d$time_mo)
    d$behav_resid <- NA_real_
    if (sum(ok) >= 10) {
      form_resid <- if (include_time) {
        behav ~ gdss_score + age + sex + time_mo
      } else {
        behav ~ gdss_score + age + sex
      }
      d$behav_resid[ok] <- resid(lm(form_resid,
                                    data = d[ok, , drop = FALSE]))
    }
    d
  }) %>%
  ungroup()

## ===================== 2) FIT MODERATION MODELS =====================
cat("--- Fitting lesion size moderation models ---\n")

all_results <- list()

for (circ in thalamic_circuits) {
  
  circ_label     <- circ$label
  circ_beta_main <- circ$beta_main
  circ_fdr_q     <- circ$fdr_q
  
  dfp <- dat_stroke_ls %>%
    filter(
      feature_type == circ$feature_type,
      feature_name == circ$feature_name,
      !is.na(value),
      !is.na(behav_resid),
      !is.na(lesion_size_z),
      !is.na(gdss_score),
      !is.na(age),
      !is.na(sex)
    ) %>%
    { if (include_time) filter(., !is.na(time_mo)) else . } %>%
    droplevels()
  
  n_obs  <- nrow(dfp)
  n_subj <- n_distinct(dfp$study_id)
  
  cat(sprintf("  %s: n_obs = %d, n_subjects = %d\n",
              circ_label, n_obs, n_subj))
  
  if (n_obs < 8 || n_subj < 5) {
    cat(sprintf("  SKIPPED: insufficient data\n"))
    next
  }
  
  # Moderation model: behav_resid * lesion_size_z
  form_mod <- if (include_time) {
    as.formula("value ~ behav_resid * lesion_size_z + gdss_score + age + sex + time_mo + (1 | study_id)")
  } else {
    as.formula("value ~ behav_resid * lesion_size_z + gdss_score + age + sex + (1 | study_id)")
  }
  
  if (nlevels(droplevels(dfp$sex)) < 2)
    form_mod <- update(form_mod, . ~ . - sex)
  
  m_mod <- try(lmer(form_mod, data = dfp, REML = FALSE, control = ctrl),
               silent = TRUE)
  
  if (inherits(m_mod, "try-error")) {
    cat(sprintf("  FAILED: %s\n",
                conditionMessage(attr(m_mod, "condition"))))
    next
  }
  
  td <- broom.mixed::tidy(m_mod, effects = "fixed") %>%
    ensure_pvalues() %>%
    mutate(
      circuit             = circ_label,
      n_obs               = n_obs,
      n_subjects          = n_subj,
      main_beta_published = circ_beta_main,
      main_fdr_published  = circ_fdr_q
    )
  
  all_results[[length(all_results) + 1]] <- td
}

results_df <- bind_rows(all_results)

## ===================== 3) SUMMARY TABLE =====================
summary_table <- results_df %>%
  filter(term %in% c("behav_resid", "lesion_size_z",
                     "behav_resid:lesion_size_z")) %>%
  select(circuit, term, estimate, std.error, statistic,
         p.value, n_obs, n_subjects,
         main_beta_published, main_fdr_published) %>%
  mutate(p_label = p_fmt(p.value)) %>%
  arrange(circuit, term)

cat("\n=== LESION SIZE MODERATION SUMMARY ===\n")
print(summary_table %>%
        select(circuit, term, estimate, p_label, n_subjects) %>%
        as.data.frame(),
      row.names = FALSE)

# Key interaction terms
interaction_terms <- summary_table %>%
  filter(term == "behav_resid:lesion_size_z")

cat("\n=== KEY: behav_resid x lesion_size_z interaction ===\n")
print(interaction_terms %>%
        select(circuit, estimate, std.error, p_label, n_subjects) %>%
        as.data.frame(),
      row.names = FALSE)

cat("\n=== INTERPRETATION ===\n")
for (i in seq_len(nrow(interaction_terms))) {
  r <- interaction_terms[i, ]
  direction <- ifelse(r$estimate > 0,
                      "larger lesions WEAKEN the fatigue-connectivity relationship",
                      "larger lesions STRENGTHEN the fatigue-connectivity relationship")
  sig <- ifelse(as.numeric(r$p.value) < 0.05, "SIGNIFICANT", "non-significant")
  cat(sprintf("  %s: %s (p = %s) — %s\n",
              r$circuit, sig, r$p_label, direction))
}
cat("\n")

## ===================== 4) SIMPLE SLOPES PLOT =====================
# Visualize the interaction at low/medium/high lesion size
# (if interaction is significant or trending)

for (circ in thalamic_circuits) {
  
  dfp <- dat_stroke_ls %>%
    filter(
      feature_type == circ$feature_type,
      feature_name == circ$feature_name,
      !is.na(value), !is.na(behav_resid),
      !is.na(lesion_size_z), !is.na(gdss_score),
      !is.na(age), !is.na(sex)
    ) %>%
    { if (include_time) filter(., !is.na(time_mo)) else . }
  
  if (nrow(dfp) < 10) next
  
  # Residualize connectivity for covariates
  form_fc_resid <- if (include_time) {
    value ~ gdss_score + age + sex + time_mo + (1 | study_id)
  } else {
    value ~ gdss_score + age + sex + (1 | study_id)
  }
  
  m_cov <- try(lmer(form_fc_resid, data = dfp, REML = FALSE,
                    control = ctrl), silent = TRUE)
  if (inherits(m_cov, "try-error")) next
  
  dfp <- dfp %>% mutate(fc_resid = resid(m_cov))
  
  # Tertile split for visualization
  dfp <- dfp %>%
    mutate(
      lesion_tertile = cut(
        lesion_size_log,
        breaks = quantile(lesion_size_log,
                          probs = c(0, 1/3, 2/3, 1),
                          na.rm = TRUE),
        labels = c("Small", "Medium", "Large"),
        include.lowest = TRUE
      )
    )
  
  # Correlation within each tertile
  cors <- dfp %>%
    filter(!is.na(lesion_tertile)) %>%
    group_by(lesion_tertile) %>%
    summarise(
      r = cor(behav_resid, fc_resid, use = "complete.obs"),
      n = n_distinct(study_id),
      .groups = "drop"
    )
  
  cat(sprintf("Simple slopes by lesion size tertile — %s:\n", circ$label))
  print(cors)
  cat("\n")
  
  p_scatter <- ggplot(
    dfp %>% filter(!is.na(lesion_tertile)),
    aes(x = behav_resid, y = fc_resid, color = lesion_tertile)
  ) +
    geom_point(alpha = 0.5, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
    scale_color_manual(
      values = c("Small" = "#185FA5", "Medium" = "#888780",
                 "Large" = "#D85A30"),
      name   = "Lesion size"
    ) +
    labs(
      title    = circ$label,
      subtitle = "Fatigue-connectivity slope by lesion size tertile",
      x        = "Vitality residual (lower = greater fatigue)",
      y        = "Connectivity residual"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(out_dir, sprintf("LESIONSIZE_%s_simple_slopes.png",
                               gsub("[^A-Za-z0-9]", "_", circ$label))),
    p_scatter, width = 6, height = 5, dpi = 300, bg = "white"
  )
  cat(sprintf("  Saved simple slopes plot for %s\n\n", circ$label))
}

## ===================== 5) SAVE OUTPUTS =====================
readr::write_csv(
  results_df,
  file.path(out_dir, "LESIONSIZE_thalamic_moderation_full.csv")
)

readr::write_csv(
  summary_table,
  file.path(out_dir, "LESIONSIZE_thalamic_moderation_summary.csv")
)

cat("\u2705 Lesion size moderation outputs saved:\n")
cat("   LESIONSIZE_thalamic_moderation_full.csv\n")
cat("   LESIONSIZE_thalamic_moderation_summary.csv\n")
cat("   LESIONSIZE_*_simple_slopes.png (one per circuit)\n")

cat("\n=== REPORTING TEMPLATE ===\n")
for (i in seq_len(nrow(interaction_terms))) {
  r <- interaction_terms[i, ]
  cat(sprintf(
    "  %s: vitality x lesion size interaction: β = %.4f, p = %s (n = %d subjects)\n",
    r$circuit, r$estimate, r$p_label, r$n_subjects
  ))
}
