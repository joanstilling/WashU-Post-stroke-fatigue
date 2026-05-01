suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lme4)
  library(broom.mixed)
})

## ---- CONFIG: EDIT THESE TWO PATHS ----
merged_path <- "/Users/.../merged_imaging_behaviour_long_ALL.csv"
out_dir     <- "/Users/..."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha_fdr     <- 0.05
include_time  <- TRUE  # derive time_mo if missing

feature_types_keep <- c(
  "WB_FC","NET_BETWEEN_FC","NET_WITHIN_FC",
  "DLPFC_TO_NET_FC","DLPFC_TO_ROI_FC",
  "SEED_EDGE_FC","SEED_TO_NET_FC",
  "CHACO_NET_MEAN","CHACO_NET_WITHIN_MEAN","CHACO_NET_BETWEEN_MEAN"
)

## ---- Helpers ----
derive_time_mo <- function(event_name) {
  e <- tolower(as.character(event_name))
  dplyr::case_when(
    e %in% c("acute_arm_1") ~ 0,
    e %in% c("3month_arm_1", "visit_1_arm_2") ~ 3,
    e %in% c("1year_arm_1", "visit_2_arm_2") ~ 12,
    TRUE ~ NA_real_
  )
}

ensure_pvalues <- function(td) {
  if (!("statistic" %in% names(td)) && ("t.value" %in% names(td))) {
    td <- td %>% rename(statistic = t.value)
  }
  if (!("p.value" %in% names(td)) && ("statistic" %in% names(td))) {
    td <- td %>% mutate(p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE))
  }
  td
}

add_model_family <- function(feature_type) {
  ft <- tolower(trimws(as.character(feature_type)))
  dplyr::case_when(
    ft %in% c("wb_fc","net_within_fc","net_between_fc") ~ "FC_NETWORK",
    ft %in% c("dlpfc_to_net_fc","dlpfc_to_roi_fc") ~ "FC_DLPFC",
    ft %in% c("seed_edge_fc","seed_to_net_fc","seed_to_roi_fc") ~ "FC_SEED",
    ft %in% c("chaco_net_mean","chaco_net_within_mean","chaco_net_between_mean") ~ "CHACO",
    stringr::str_detect(ft, "chaco") ~ "CHACO",
    TRUE ~ "UNMAPPED"
  )
}

safe_file <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

############################################################
## 1) LOAD MERGED CSV
############################################################
dat <- readr::read_csv(merged_path, show_col_types = FALSE)

cat("\nLoaded:", merged_path, "\n")
cat("Rows:", nrow(dat), " Unique subjects:", dplyr::n_distinct(dat$study_id), "\n")

# ✅ ADD gdss_score to required columns
required_cols <- c("study_id","event_name","feature_type","feature_name","value",
                   "behav","gdss_score","age","sex","group")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Merged CSV missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("\nUnique feature_type values in your CSV (this is the truth):\n")
print(sort(unique(dat$feature_type)))

############################################################
## 2) CLEAN + FILTER + COERCE TYPES + RECODE GROUP/SEX
############################################################
dat_model <- dat %>%
  mutate(
    study_id      = as.character(study_id),
    event_name    = as.character(event_name),
    feature_type  = as.character(feature_type),
    feature_name  = as.character(feature_name),
    value = suppressWarnings(as.numeric(value)),
    behav = suppressWarnings(as.numeric(behav)),
    gdss_score = suppressWarnings(as.numeric(gdss_score)),
    age   = suppressWarnings(as.numeric(age))
  )

cat("\nRows with 'chaco' in feature_type (before keep filter):\n")
print(dat_model %>% filter(str_detect(tolower(feature_type), "chaco")) %>% count(feature_type, sort = TRUE))

dat_model <- dat_model %>% filter(feature_type %in% feature_types_keep)

cat("\nAfter feature_types_keep filter:\n")
print(table(dat_model$feature_type))

if (include_time) {
  if (!("time_mo" %in% names(dat_model))) {
    dat_model <- dat_model %>% mutate(time_mo = derive_time_mo(event_name))
    cat("\nDerived time_mo from event_name.\n")
  } else {
    dat_model <- dat_model %>% mutate(time_mo = suppressWarnings(as.numeric(time_mo)))
    cat("\nUsing time_mo from CSV.\n")
  }
}

# ✅ Your coding: group 0=stroke, 1=control; sex 0=female, 1=male
dat_model <- dat_model %>%
  mutate(
    group = case_when(
      group %in% c("0", 0) ~ "stroke",
      group %in% c("1", 1) ~ "control",
      TRUE ~ as.character(group)
    ),
    sex = case_when(
      sex %in% c("0", 0) ~ "female",
      sex %in% c("1", 1) ~ "male",
      TRUE ~ as.character(sex)
    ),
    group = factor(group),
    sex   = factor(sex),
    model_family = add_model_family(feature_type)
  )

if ("control" %in% levels(dat_model$group)) dat_model$group <- relevel(dat_model$group, ref = "control")
if ("female"  %in% levels(dat_model$sex))   dat_model$sex   <- relevel(dat_model$sex,   ref = "female")

cat("\nModel_family counts (after keep filter):\n")
print(table(dat_model$model_family, useNA = "ifany"))
############################################################
## DIAGNOSTICS: correlation + VIF (collinearity)
############################################################

# collapse to unique subject x event (or time_mo) BEFORE correlating
cor_visit <- dat_model %>%
  distinct(study_id, event_name, group, behav, gdss_score) %>%
  filter(is.finite(behav), is.finite(gdss_score)) %>%
  group_by(group) %>%
  summarise(
    n_subject_visits = n(),
    pearson  = cor(behav, gdss_score, use="complete.obs"),
    spearman = cor(behav, gdss_score, use="complete.obs", method="spearman"),
    .groups="drop"
  )

print(cor_visit)

# optional: subject-level (average across visits)
cor_subject <- dat_model %>%
  distinct(study_id, group, event_name, behav, gdss_score) %>%
  group_by(study_id, group) %>%
  summarise(
    behav_mean = mean(behav, na.rm=TRUE),
    gdss_mean  = mean(gdss_score, na.rm=TRUE),
    .groups="drop"
  ) %>%
  group_by(group) %>%
  summarise(
    n_subjects = n(),
    pearson  = cor(behav_mean, gdss_mean, use="complete.obs"),
    spearman = cor(behav_mean, gdss_mean, use="complete.obs", method="spearman"),
    .groups="drop"
  )

print(cor_subject)

vif_like <- function(fit) {
  X <- model.matrix(fit)
  cn <- colnames(X)
  keep <- cn != "(Intercept)"
  X <- X[, keep, drop=FALSE]
  out <- sapply(seq_len(ncol(X)), function(j) {
    y <- X[, j]; Xo <- X[, -j, drop=FALSE]
    r2 <- summary(lm(y ~ Xo))$r.squared
    1/(1-r2)
  })
  data.frame(term = cn[keep], vif = out, row.names=NULL)
}

design <- dat_model %>%
  distinct(study_id, event_name, group, age, sex, time_mo, behav, gdss_score) %>%
  filter(is.finite(behav), is.finite(gdss_score), is.finite(age), !is.na(sex), !is.na(group),
         if (include_time) is.finite(time_mo) else TRUE)

fit_design <- if (include_time) {
  lm(behav ~ gdss_score + group + age + sex + time_mo, data = design)
} else {
  lm(behav ~ gdss_score + group + age + sex, data = design)
}
print(vif_like(fit_design))

############################################################
# Borrow behaviour + GDSS across timepoints ONLY for CHACO
# (uses study_id now)
############################################################

## Ensure time_mo exists (derive if missing)
if (!("time_mo" %in% names(dat_model))) {
  dat_model <- dat_model %>% mutate(time_mo = derive_time_mo(event_name))
} else {
  dat_model <- dat_model %>% mutate(time_mo = suppressWarnings(as.numeric(time_mo)))
}

# Create lookups (subject x time_mo) from non-missing values
behav_lookup <- dat_model %>%
  filter(!is.na(behav), !is.na(time_mo)) %>%
  distinct(study_id, time_mo, behav)

gdss_lookup <- dat_model %>%
  filter(!is.na(gdss_score), !is.na(time_mo)) %>%
  distinct(study_id, time_mo, gdss_score)

# Fill missing behav/GDSS ONLY for CHACO rows with nearest available value
dat_model <- dat_model %>%
  mutate(is_chaco = str_detect(tolower(feature_type), "chaco")) %>%
  group_by(study_id) %>%
  group_modify(function(df, key) {
    
    if (!any(df$is_chaco, na.rm = TRUE)) return(df)
    
    subj_beh <- behav_lookup %>% filter(study_id == key$study_id)
    subj_gds <- gdss_lookup  %>% filter(study_id == key$study_id)
    
    idx_fill <- which(df$is_chaco & !is.na(df$time_mo) &
                        (is.na(df$behav) | is.na(df$gdss_score)))
    if (length(idx_fill) == 0) return(df)
    
    for (ii in idx_fill) {
      t0 <- df$time_mo[ii]
      
      # --- fill behav if missing ---
      if (is.na(df$behav[ii]) && nrow(subj_beh) > 0) {
        cand <- subj_beh %>% mutate(dt = abs(time_mo - t0), tie = time_mo) %>% arrange(dt, tie)
        if (nrow(cand) > 0) df$behav[ii] <- cand$behav[1]
      }
      
      # --- fill gdss_score if missing ---
      if (is.na(df$gdss_score[ii]) && nrow(subj_gds) > 0) {
        cand <- subj_gds %>% mutate(dt = abs(time_mo - t0), tie = time_mo) %>% arrange(dt, tie)
        if (nrow(cand) > 0) df$gdss_score[ii] <- cand$gdss_score[1]
      }
    }
    df
  }) %>%
  ungroup() %>%
  select(-is_chaco)

audit_chaco_fill <- dat_model %>%
  mutate(is_chaco = str_detect(tolower(feature_type), "chaco")) %>%
  summarise(
    chaco_rows = sum(is_chaco),
    chaco_missing_behav_after = sum(is_chaco & is.na(behav)),
    chaco_missing_gdss_after  = sum(is_chaco & is.na(gdss_score)),
    chaco_complete_after = sum(is_chaco & !is.na(value) & !is.na(behav) & !is.na(gdss_score) &
                                 !is.na(age) & !is.na(sex) & !is.na(group) & !is.na(time_mo))
  )
print(audit_chaco_fill)

############################################################
## 2b) CHECK: Are CHACO rows being dropped by missingness?
############################################################
miss_by_family <- dat_model %>%
  mutate(
    complete = !is.na(value) & !is.na(behav) & !is.na(gdss_score) &
      !is.na(age) & !is.na(sex) & !is.na(group) &
      (if (include_time) !is.na(time_mo) else TRUE)
  ) %>%
  group_by(model_family) %>%
  summarise(
    n_rows = n(),
    n_complete = sum(complete),
    n_subjects = n_distinct(study_id),
    n_subjects_complete = n_distinct(study_id[complete]),
    .groups = "drop"
  )

cat("\nCompleteness by family (now requires GDSS too for 'complete'):\n")
print(miss_by_family)

############################################################
# ChaCo moderation mixed-effects (STROKE ONLY)
# Now adds GDSS as covariate (z-scored within stroke)
############################################################

chaco_pick_feature_type <- c("CHACO_NET_MEAN", "CHACO_NET_WITHIN_MEAN", "CHACO_NET_BETWEEN_MEAN")
chaco_pick_feature_name <- NULL
moderation_families <- c("FC_NETWORK", "FC_DLPFC", "FC_SEED")
try_random_slope <- TRUE

chaco_subject <- dat_model %>%
  filter(group == "stroke") %>%
  filter(feature_type %in% chaco_pick_feature_type) %>%
  filter(!is.na(value)) %>%
  {
    if (!is.null(chaco_pick_feature_name)) {
      dplyr::filter(., feature_name == chaco_pick_feature_name)
    } else .
  } %>%
  arrange(study_id, feature_type, feature_name, time_mo) %>%
  group_by(study_id) %>%
  summarise(
    chaco_value = mean(value, na.rm = TRUE),
    n_chaco_rows = n(),
    .groups = "drop"
  )

cat("\n[ChaCo moderation] Stroke subjects with ChaCo:", nrow(chaco_subject), "\n")

dat_model2 <- dat_model %>%
  left_join(chaco_subject, by = "study_id")

dat_stroke_mod <- dat_model2 %>%
  filter(group == "stroke") %>%
  mutate(
    behav_z = as.numeric(scale(behav)),
    chaco_z = as.numeric(scale(chaco_value)),
    gdss_z  = as.numeric(scale(gdss_score))
  ) %>%
  filter(
    !is.na(value),
    !is.na(behav_z),
    !is.na(chaco_z),
    !is.na(gdss_z),
    !is.na(age),
    !is.na(sex),
    !is.na(study_id),
    if (include_time) !is.na(time_mo) else TRUE
  ) %>%
  mutate(model_family = add_model_family(feature_type)) %>%
  filter(model_family %in% moderation_families) %>%
  droplevels()

cat("\n[ChaCo moderation] Stroke-only moderation dataset:\n")
cat("  Rows:", nrow(dat_stroke_mod),
    " Subjects:", dplyr::n_distinct(dat_stroke_mod$study_id),
    " Features:", nrow(distinct(dat_stroke_mod, feature_type, feature_name)), "\n")

form_int_only <- if (include_time) {
  as.formula("value ~ behav_z * chaco_z + gdss_z + age + sex + time_mo + (1 | study_id)")
} else {
  as.formula("value ~ behav_z * chaco_z + gdss_z + age + sex + (1 | study_id)")
}

form_rand_slope <- if (include_time) {
  as.formula("value ~ behav_z * chaco_z + gdss_z + age + sex + time_mo + (1 + behav_z | study_id)")
} else {
  as.formula("value ~ behav_z * chaco_z + gdss_z + age + sex + (1 + behav_z | study_id)")
}

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

results_chaco_mod <- dat_stroke_mod %>%
  group_by(feature_type, feature_name) %>%
  group_modify(function(df, key) {
    
    n_subj <- n_distinct(df$study_id)
    if (nrow(df) < 12 || n_subj < 8) return(tibble())
    
    df2 <- df %>% droplevels()
    f_int <- form_int_only
    f_slp <- form_rand_slope
    if (nlevels(df2$sex) < 2) {
      f_int <- update(f_int, . ~ . - sex)
      f_slp <- update(f_slp, . ~ . - sex)
    }
    
    m <- NULL
    used_random_slope <- FALSE
    if (isTRUE(try_random_slope)) {
      m_try <- try(lmer(f_slp, data = df2, REML = FALSE, control = ctrl), silent = TRUE)
      if (!inherits(m_try, "try-error")) { m <- m_try; used_random_slope <- TRUE }
    }
    if (is.null(m)) {
      m_try2 <- try(lmer(f_int, data = df2, REML = FALSE, control = ctrl), silent = TRUE)
      if (!inherits(m_try2, "try-error")) { m <- m_try2; used_random_slope <- FALSE }
    }
    if (is.null(m)) return(tibble())
    
    td <- broom.mixed::tidy(m, effects = "fixed") %>% ensure_pvalues()
    
    td %>%
      transmute(
        term, estimate, std.error,
        statistic = if ("statistic" %in% names(.)) statistic else NA_real_,
        p.value   = if ("p.value" %in% names(.)) p.value else NA_real_,
        n_obs = nrow(df2),
        n_subjects = n_distinct(df2$study_id),
        model_family = add_model_family(key$feature_type),
        used_random_slope = used_random_slope,
        effect_type = case_when(
          term == "behav_z" ~ "vitality_main",
          term == "gdss_z"  ~ "gdss_main",
          term == "chaco_z" ~ "chaco_main",
          term %in% c("behav_z:chaco_z", "chaco_z:behav_z") ~ "vitality_by_chaco",
          TRUE ~ "other"
        )
      )
  }) %>%
  ungroup()

readr::write_csv(
  results_chaco_mod,
  file.path(out_dir, "STROKEONLY_ChacoModeration_all_feature_models_mixed.csv")
)

mod_terms <- results_chaco_mod %>%
  filter(effect_type == "vitality_by_chaco", !is.na(p.value)) %>%
  group_by(model_family) %>%
  mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fdr)

readr::write_csv(
  mod_terms,
  file.path(out_dir, "STROKEONLY_ChacoModeration_vitality_by_chaco_with_FDR.csv")
)
############################################################
## Residualize vitality on depression (within feature):
## behav_resid = behav with GDSS-related variance removed
############################################################
dat_model <- dat_model %>%
  group_by(feature_type, feature_name) %>%
  group_modify(~{
    d <- .x
    ok <- is.finite(d$behav) & is.finite(d$gdss_score)
    d$behav_resid <- NA_real_
    if (sum(ok) >= 10) {
      d$behav_resid[ok] <- resid(lm(behav ~ gdss_score, data = d[ok, , drop = FALSE]))
    }
    d
  }) %>%
  ungroup()

############################################################
## 3) FIT PER-FEATURE MIXED MODELS (ALL SUBJECTS)
## NOW: value ~ behav * group + gdss_score + covariates + (1|study_id)
############################################################
base_form <- if (include_time) {
  as.formula("value ~ behav * group + gdss_score + age + sex + time_mo + (1 | study_id)")
} else {
  as.formula("value ~ behav * group + gdss_score + age + sex + (1 | study_id)")
}

# --- CHACO survival now checks both behav + gdss_score ---
chaco_feats <- dat_model %>% filter(model_family == "CHACO") %>% distinct(feature_type, feature_name)
cat("Unique CHACO features:", nrow(chaco_feats), "\n")

chaco_survival <- dat_model %>%
  filter(model_family == "CHACO") %>%
  group_by(feature_type, feature_name) %>%
  summarise(
    n_rows_total = n(),
    n_rows_cc = sum(
      !is.na(value) & !is.na(behav) & !is.na(gdss_score) & !is.na(age) & !is.na(sex) & !is.na(group) &
        !is.na(study_id) & (if (include_time) !is.na(time_mo) else TRUE)
    ),
    n_subj_cc = n_distinct(study_id[
      !is.na(value) & !is.na(behav) & !is.na(gdss_score) & !is.na(age) & !is.na(sex) & !is.na(group) &
        !is.na(study_id) & (if (include_time) !is.na(time_mo) else TRUE)
    ]),
    n_groups_cc = n_distinct(group[
      !is.na(value) & !is.na(behav) & !is.na(gdss_score) & !is.na(age) & !is.na(sex) & !is.na(group) &
        !is.na(study_id) & (if (include_time) !is.na(time_mo) else TRUE)
    ]),
    .groups = "drop"
  ) %>%
  arrange(n_rows_cc, n_subj_cc)

print(head(chaco_survival, 30))

results <- dat_model %>%
  group_by(feature_type, feature_name) %>%
  group_modify(function(df, key) {
    
    df2 <- df %>%
      filter(
        !is.na(value),
        !is.na(behav),
        !is.na(gdss_score),
        !is.na(age),
        !is.na(sex),
        !is.na(group),
        !is.na(study_id),
        if (include_time) !is.na(time_mo) else TRUE
      ) %>%
      droplevels()
    
    n_subj <- n_distinct(df2$study_id)
    if (nrow(df2) < 8 || n_subj < 3) return(tibble())
    
    form_use <- base_form
    
    # If only one group, drop group + interaction but keep gdss_score
    if (nlevels(df2$group) < 2) {
      form_use <- if (include_time) {
        as.formula("value ~ behav + gdss_score + age + sex + time_mo + (1 | study_id)")
      } else {
        as.formula("value ~ behav + gdss_score + age + sex + (1 | study_id)")
      }
    }
    if (nlevels(df2$sex) < 2) {
      form_use <- update(form_use, . ~ . - sex)
    }
    
    m <- try(lmer(form_use, data = df2, REML = FALSE), silent = TRUE)
    if (inherits(m, "try-error")) return(tibble())
    
    td <- broom.mixed::tidy(m, effects = "fixed") %>% ensure_pvalues()
    
    td %>%
      transmute(
        term,
        estimate,
        std.error,
        statistic = if ("statistic" %in% names(.)) statistic else NA_real_,
        p.value   = if ("p.value" %in% names(.)) p.value else NA_real_,
        n_obs = nrow(df2),
        n_subjects = n_subj
      )
    
  }) %>%
  ungroup()

cat("\nFinished models. Feature types in results:\n")
print(table(results$feature_type))
############################################################
## SECOND RUN: Residualized vitality model (reviewer-proof)
## value ~ behav_resid * group + gdss_score + covariates + (1|study_id)
############################################################
base_form_resid <- if (include_time) {
  as.formula("value ~ behav_resid * group + gdss_score + age + sex + time_mo + (1 | study_id)")
} else {
  as.formula("value ~ behav_resid * group + gdss_score + age + sex + (1 | study_id)")
}

results_resid <- dat_model %>%
  group_by(feature_type, feature_name) %>%
  group_modify(function(df, key) {
    
    df2 <- df %>%
      filter(
        !is.na(value),
        !is.na(behav_resid),
        !is.na(gdss_score),
        !is.na(age),
        !is.na(sex),
        !is.na(group),
        !is.na(study_id),
        if (include_time) !is.na(time_mo) else TRUE
      ) %>%
      droplevels()
    
    n_subj <- n_distinct(df2$study_id)
    if (nrow(df2) < 8 || n_subj < 3) return(tibble())
    
    form_use <- base_form_resid
    
    # If only one group, drop group + interaction (keep gdss_score)
    if (nlevels(df2$group) < 2) {
      form_use <- if (include_time) {
        as.formula("value ~ behav_resid + gdss_score + age + sex + time_mo + (1 | study_id)")
      } else {
        as.formula("value ~ behav_resid + gdss_score + age + sex + (1 | study_id)")
      }
    }
    
    # If sex has 1 level in this subset, drop sex
    if (nlevels(df2$sex) < 2) {
      form_use <- update(form_use, . ~ . - sex)
    }
    
    m <- try(lmer(form_use, data = df2, REML = FALSE), silent = TRUE)
    if (inherits(m, "try-error")) return(tibble())
    
    td <- broom.mixed::tidy(m, effects = "fixed")
    td <- ensure_pvalues(td)
    
    td %>%
      transmute(
        term,
        estimate,
        std.error,
        statistic = if ("statistic" %in% names(.)) statistic else NA_real_,
        p.value   = if ("p.value" %in% names(.)) p.value else NA_real_,
        n_obs = nrow(df2),
        n_subjects = n_subj
      )
  }) %>%
  ungroup()

# Save residualized run
readr::write_csv(results_resid, file.path(out_dir, "ALL_all_feature_models_mixed_BEHAV_RESID.csv"))

############################################################
## FDR for residualized vitality terms
############################################################
results_resid_tagged <- results_resid %>%
  mutate(
    model_family = add_model_family(feature_type),
    effect_type = case_when(
      term == "behav_resid" ~ "vitality_resid_main",
      str_starts(term, "behav_resid:group") ~ "vitality_resid_by_group",
      TRUE ~ "other"
    )
  )

resid_terms_fdr <- results_resid_tagged %>%
  filter(effect_type != "other", !is.na(p.value)) %>%
  group_by(model_family, effect_type) %>%
  mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_resid_main <- resid_terms_fdr %>% filter(effect_type == "vitality_resid_main", p_fdr < alpha_fdr)
sig_resid_int  <- resid_terms_fdr %>% filter(effect_type == "vitality_resid_by_group", p_fdr < alpha_fdr)

readr::write_csv(resid_terms_fdr, file.path(out_dir, "ALL_vitality_resid_terms_with_FDR_mixed.csv"))
readr::write_csv(sig_resid_main,  file.path(out_dir, "ALL_significant_vitality_resid_main_FDR.csv"))
readr::write_csv(sig_resid_int,   file.path(out_dir, "ALL_significant_vitality_resid_by_group_FDR.csv"))

results_tagged <- results %>%
  mutate(
    model_family = add_model_family(feature_type),
    effect_type = case_when(
      term == "behav" ~ "vitality_main",
      str_starts(term, "behav:group") ~ "vitality_by_group",
      term == "gdss_score" ~ "gdss_main",
      TRUE ~ "other"
    )
  )

write_csv(results, file.path(out_dir, "ALL_all_feature_models_mixed.csv"))

############################################################
## 4) FDR WITHIN FAMILY + SAVE PER FAMILY (now includes gdss_main too)
############################################################
behav_terms_fdr <- results_tagged %>%
  filter(effect_type != "other", !is.na(p.value)) %>%
  group_by(model_family, effect_type) %>%
  mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_vit_main <- behav_terms_fdr %>% filter(effect_type == "vitality_main", p_fdr < alpha_fdr)
sig_vit_int  <- behav_terms_fdr %>% filter(effect_type == "vitality_by_group", p_fdr < alpha_fdr)
sig_gdss     <- behav_terms_fdr %>% filter(effect_type == "gdss_main", p_fdr < alpha_fdr)

write_csv(behav_terms_fdr, file.path(out_dir, "ALL_terms_with_FDR_mixed.csv"))
write_csv(sig_vit_main,    file.path(out_dir, "ALL_significant_vitality_main_effects_mixed.csv"))
write_csv(sig_vit_int,     file.path(out_dir, "ALL_significant_vitality_by_group_effects_mixed.csv"))
write_csv(sig_gdss,        file.path(out_dir, "ALL_significant_gdss_main_effects_mixed.csv"))

families_present <- sort(unique(results_tagged$model_family))
cat("\nSaving families:\n"); print(families_present)

for (fam in families_present) {
  fam_tag <- safe_file(fam)
  
  res_all  <- results_tagged  %>% filter(model_family == fam)
  res_fdr  <- behav_terms_fdr %>% filter(model_family == fam)
  res_vit  <- sig_vit_main    %>% filter(model_family == fam)
  res_int  <- sig_vit_int     %>% filter(model_family == fam)
  res_gds  <- sig_gdss        %>% filter(model_family == fam)
  
  write_csv(res_all, file.path(out_dir, paste0(fam_tag, "_all_feature_models_mixed.csv")))
  write_csv(res_fdr, file.path(out_dir, paste0(fam_tag, "_terms_with_FDR_mixed.csv")))
  write_csv(res_vit, file.path(out_dir, paste0(fam_tag, "_sig_vitality_main_FDR.csv")))
  write_csv(res_int, file.path(out_dir, paste0(fam_tag, "_sig_vitality_by_group_FDR.csv")))
  write_csv(res_gds, file.path(out_dir, paste0(fam_tag, "_sig_gdss_main_FDR.csv")))
  
  cat("Saved ", fam, ": ", nrow(res_all), " rows\n", sep = "")
}

cat("\n✅ Done (with GDSS). Outputs saved to:\n", out_dir, "\n")
