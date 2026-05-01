
# =============================================================================
#  SUPPLEMENTARY REGRESSION TABLES
#  Table S1: Depression-adjusted primary model
#  Table S2: Residualized vitality sensitivity analysis
#  Table S3: ChaCo structural disconnection moderation
#  Table S4: FC node strength group differences (stroke vs control, GDSS-adj)
#  Table S5: ALFF/fALFF group differences (stroke vs control, GDSS-adj)
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(openxlsx)
})

# --- paths -------------------------------------------------------------------
primary_path      <- "/Users/.../ALL_terms_with_FDR_mixed.csv"
primary_alff_path <- "/Users/.../ALL_ALFF_fALFF_terms_with_FDR_mixed.csv"
resid_path        <- "/Users/.../ALL_vitality_resid_terms_with_FDR_mixed.csv"
resid_alff_path   <- "/Users/.../ALL_ALFF_fALFF_terms_with_FDR_mixed_BEHAV_RESID.csv"
chaco_fc_path     <- "/Users/.../STROKEONLY_ChacoModeration_vitality_by_chaco_with_FDR.csv"
chaco_alff_path   <- "/Users/.../STROKEONLY_ChacoModeration_vitality_by_chaco_with_FDR_ALFF_fALFF.csv"
alff_falff_path   <- "/Users/.../ALLTIME_ALL_group_tvals_GDSSadj_combined.csv"
out_path          <- "/Users/.../supplementary_tables.xlsx"

# Filter threshold for regression tables
p_threshold <- 0.05

# --- helpers -----------------------------------------------------------------
fmt_beta  <- function(x) sprintf("%.4f", x)
fmt_se    <- function(x) sprintf("%.4f", x)
fmt_t     <- function(x) sprintf("%.3f", x)
fmt_p     <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
fmt_q     <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))

clean_feature <- function(x) {
  x %>%
    str_remove("^SEED:") %>%
    str_remove("^NS_") %>%
    str_replace_all("->", " \u2192 ") %>%
    str_replace_all("<->", " \u2194 ") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("THalamus", "Thalamus") %>%
    str_replace_all("CaudalACC", "Caudal ACC") %>%
    str_replace_all("ctx-lh-", "LH ") %>%
    str_replace_all("ctx-rh-", "RH ")
}

clean_term <- function(x) {
  case_when(
    x == "behav"                    ~ "Vitality",
    x == "behav:groupstroke"        ~ "Vitality x Group (stroke)",
    x == "behav_resid"              ~ "Vitality (residualized)",
    x == "behav_resid:groupstroke"  ~ "Vitality (residualized) x Group (stroke)",
    x == "gdss_score"               ~ "GDSS (depressive symptoms)",
    TRUE ~ x
  )
}

family_order <- c(
  "FC_NETWORK",
  "FC_SEED",
  "FC_DLPFC",
  "ALFF_NETWORK",
  "ALFF_ROI",
  "FALFF_NETWORK",
  "FALFF_ROI",
  "CHACO"
)

family_label <- function(x) {
  case_when(
    x == "FC_NETWORK"    ~ "Between-network functional connectivity",
    x == "FC_DLPFC"      ~ "DLPFC ROI-based connectivity",
    x == "FC_SEED"       ~ "Thalamic ROI-based connectivity",
    x == "ALFF_ROI"      ~ "ALFF (regional amplitude)",
    x == "FALFF_ROI"     ~ "fALFF (fractional amplitude)",
    x == "ALFF_NETWORK"  ~ "ALFF (network)",
    x == "FALFF_NETWORK" ~ "fALFF (network)",
    x == "CHACO"         ~ "Structural disconnection (ChaCo)",
    TRUE ~ x
  )
}

sig_star <- function(p_fdr, p_raw) {
  case_when(
    p_fdr < 0.05  ~ "*",
    p_raw < 0.05  ~ "+",
    TRUE          ~ ""
  )
}

# --- styles ------------------------------------------------------------------
make_styles <- function() {
  list(
    title = createStyle(
      fontName = "Arial", fontSize = 11, fontColour = "#2C2C2A",
      textDecoration = "bold", wrapText = FALSE
    ),
    header = createStyle(
      fontName = "Arial", fontSize = 10, fontColour = "white",
      fgFill = "#185FA5", textDecoration = "bold",
      halign = "center", valign = "center",
      border = "Bottom", borderColour = "#0C447C", borderStyle = "medium",
      wrapText = TRUE
    ),
    family = createStyle(
      fontName = "Arial", fontSize = 10, fontColour = "#2C2C2A",
      textDecoration = "bold", fgFill = "#E6F1FB",
      border = "TopBottom", borderColour = "#B5D4F4", borderStyle = "thin"
    ),
    data = createStyle(
      fontName = "Arial", fontSize = 9, fontColour = "#2C2C2A",
      halign = "left", valign = "center"
    ),
    data_center = createStyle(
      fontName = "Arial", fontSize = 9, fontColour = "#2C2C2A",
      halign = "center", valign = "center"
    ),
    sig_fdr = createStyle(
      fontName = "Arial", fontSize = 9, fontColour = "#185FA5",
      textDecoration = "bold", halign = "center"
    ),
    sig_unc = createStyle(
      fontName = "Arial", fontSize = 9, fontColour = "#888780",
      halign = "center"
    ),
    note = createStyle(
      fontName = "Arial", fontSize = 9, fontColour = "#5F5E5A",
      textDecoration = "italic"
    )
  )
}

# =============================================================================
#  FUNCTION: write regression table sheet (S1/S2/S3)
# =============================================================================
write_table_sheet <- function(wb, sheet_name, tbl, table_title, table_note, sty) {
  
  cat(sprintf("Writing sheet: %s (%d rows)\n", sheet_name, nrow(tbl)))
  addWorksheet(wb, sheet_name, gridLines = FALSE)
  
  writeData(wb, sheet_name, table_title, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name, sty$title, rows = 1, cols = 1)
  
  col_names <- c("Model family", "Feature", "Term",
                 "\u03b2", "SE", "t", "p", "FDR q", "", "n obs", "n subj")
  writeData(wb, sheet_name, as.data.frame(t(col_names)),
            startRow = 2, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, sty$header,
           rows = 2, cols = 1:11, gridExpand = TRUE)
  
  current_row <- 3
  prev_family <- ""
  
  for (i in seq_len(nrow(tbl))) {
    row_sig     <- tbl$Sig[i]
    row_family  <- tbl$Family[i]
    row_feature <- tbl$Feature[i]
    row_term    <- tbl$Term[i]
    row_beta    <- tbl$Beta[i]
    row_se      <- tbl$SE[i]
    row_t       <- tbl$t[i]
    row_p       <- tbl$p[i]
    row_q       <- tbl$`FDR q`[i]
    row_nobs    <- tbl$n_obs[i]
    row_nsubj   <- tbl$n_subj[i]
    
    if (row_family != prev_family) {
      writeData(wb, sheet_name, row_family,
                startRow = current_row, startCol = 1)
      mergeCells(wb, sheet_name, rows = current_row, cols = 1:11)
      addStyle(wb, sheet_name, sty$family,
               rows = current_row, cols = 1:11, gridExpand = TRUE)
      current_row <- current_row + 1
      prev_family <- row_family
    }
    
    writeData(wb, sheet_name,
              data.frame(V1="", V2=row_feature, V3=row_term,
                         V4=row_beta, V5=row_se, V6=row_t,
                         V7=row_p, V8=row_q, V9=row_sig,
                         V10=row_nobs, V11=row_nsubj),
              startRow = current_row, startCol = 1, colNames = FALSE)
    
    addStyle(wb, sheet_name, sty$data,
             rows = current_row, cols = 1:3, gridExpand = TRUE)
    addStyle(wb, sheet_name, sty$data_center,
             rows = current_row, cols = 4:8, gridExpand = TRUE)
    addStyle(wb, sheet_name, sty$data_center,
             rows = current_row, cols = 10:11, gridExpand = TRUE)
    
    if (row_sig == "*") {
      addStyle(wb, sheet_name, sty$sig_fdr, rows = current_row, cols = 9)
    } else if (row_sig == "+") {
      addStyle(wb, sheet_name, sty$sig_unc, rows = current_row, cols = 9)
    } else {
      addStyle(wb, sheet_name, sty$data_center, rows = current_row, cols = 9)
    }
    
    if (i %% 2 == 0) {
      addStyle(wb, sheet_name, createStyle(fgFill = "#F1EFE8"),
               rows = current_row, cols = 1:10,
               gridExpand = TRUE, stack = TRUE)
    }
    
    current_row <- current_row + 1
  }
  
  note_row <- current_row + 1
  writeData(wb, sheet_name, table_note, startRow = note_row, startCol = 1)
  mergeCells(wb, sheet_name, rows = note_row, cols = 1:10)
  addStyle(wb, sheet_name, sty$note, rows = note_row, cols = 1)
  
  setColWidths(wb, sheet_name, cols = 1:11,
               widths = c(28, 28, 30, 8, 8, 8, 8, 8, 4, 7, 7))
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 1)
}

# =============================================================================
#  FUNCTION: write group difference sheet (S4/S5)
# =============================================================================
write_group_diff_sheet <- function(wb, sheet_name, tbl, table_title,
                                   table_note, sty) {
  
  cat(sprintf("Writing sheet: %s (%d rows)\n", sheet_name, nrow(tbl)))
  addWorksheet(wb, sheet_name, gridLines = FALSE)
  
  writeData(wb, sheet_name, table_title, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name, sty$title, rows = 1, cols = 1)
  
  col_names <- c("Metric", "Region",
                 "t", "p", "FDR q", "", "Est.", "SE", "n obs", "n subj")
  writeData(wb, sheet_name, as.data.frame(t(col_names)),
            startRow = 2, startCol = 1, colNames = FALSE)
  addStyle(wb, sheet_name, sty$header,
           rows = 2, cols = 1:11, gridExpand = TRUE)
  
  current_row <- 3
  prev_metric <- ""
  
  for (i in seq_len(nrow(tbl))) {
    row_metric <- tbl$Metric[i]
    row_region <- tbl$Region[i]
    row_t      <- tbl$t[i]
    row_p      <- tbl$p[i]
    row_q      <- tbl$`FDR q`[i]
    row_sig    <- tbl$Sig[i]
    row_est    <- tbl$Est[i]
    row_se     <- tbl$SE[i]
    row_nobs   <- tbl$n_obs[i]
    row_nsubj  <- tbl$n_subj[i]
    
    # Metric subheader
    if (row_metric != prev_metric) {
      writeData(wb, sheet_name, row_metric,
                startRow = current_row, startCol = 1)
      mergeCells(wb, sheet_name, rows = current_row, cols = 1:11)
      addStyle(wb, sheet_name, sty$family,
               rows = current_row, cols = 1:11, gridExpand = TRUE)
      current_row <- current_row + 1
      prev_metric <- row_metric
    }
    
    writeData(wb, sheet_name,
              data.frame(V1="", V2=row_region,
                         V3=row_t, V4=row_p, V5=row_q, V6=row_sig,
                         V7=row_est, V8=row_se,
                         V9=row_nobs, V10=row_nsubj),
              startRow = current_row, startCol = 1, colNames = FALSE)
    
    addStyle(wb, sheet_name, sty$data,
             rows = current_row, cols = 1:2, gridExpand = TRUE)
    addStyle(wb, sheet_name, sty$data_center,
             rows = current_row, cols = 3:7, gridExpand = TRUE)
    addStyle(wb, sheet_name, sty$data_center,
             rows = current_row, cols = 9:10, gridExpand = TRUE)
    
    if (row_sig == "*") {
      addStyle(wb, sheet_name, sty$sig_fdr, rows = current_row, cols = 6)
    } else if (row_sig == "+") {
      addStyle(wb, sheet_name, sty$sig_unc, rows = current_row, cols = 6)
    } else {
      addStyle(wb, sheet_name, sty$data_center, rows = current_row, cols = 6)
    }
    
    if (i %% 2 == 0) {
      addStyle(wb, sheet_name, createStyle(fgFill = "#F1EFE8"),
               rows = current_row, cols = 1:10,
               gridExpand = TRUE, stack = TRUE)
    }
    
    current_row <- current_row + 1
  }
  
  note_row <- current_row + 1
  writeData(wb, sheet_name, table_note, startRow = note_row, startCol = 1)
  mergeCells(wb, sheet_name, rows = note_row, cols = 1:10)
  addStyle(wb, sheet_name, sty$note, rows = note_row, cols = 1)
  
  setColWidths(wb, sheet_name, cols = 1:10,
               widths = c(22, 30, 8, 8, 8, 4, 8, 8, 7, 7))
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 1)
}

# =============================================================================
#  PROCESS TABLES S1-S3
# =============================================================================
df1      <- read_csv(primary_path,      show_col_types = FALSE)
df1_alff <- read_csv(primary_alff_path, show_col_types = FALSE)
shared_cols <- intersect(names(df1), names(df1_alff))
df1_combined <- bind_rows(df1 %>% select(all_of(shared_cols)),
                          df1_alff %>% select(all_of(shared_cols)))
cat("=== PRIMARY MODEL: families present after merge ===\n")
print(table(df1_combined$model_family))

tbl1 <- df1_combined %>%
  filter(term %in% c("behav", "behav:groupstroke")) %>%
  filter(p.value < p_threshold) %>%
  mutate(family_order_idx = match(model_family, family_order)) %>%
  arrange(family_order_idx, feature_type, p.value) %>%
  mutate(Family  = family_label(model_family),
         Feature = clean_feature(feature_name),
         Term    = clean_term(term),
         Beta    = fmt_beta(estimate),
         SE      = fmt_se(std.error),
         t       = fmt_t(statistic),
         p       = fmt_p(p.value),
         `FDR q` = fmt_q(p_fdr),
         Sig     = sig_star(p_fdr, p.value),
         n_obs   = n_obs,
         n_subj  = n_subjects) %>%
  select(Family, Feature, Term, Beta, SE, t, p, `FDR q`, Sig, n_obs, n_subj)

df2      <- read_csv(resid_path,      show_col_types = FALSE)
df2_alff <- read_csv(resid_alff_path, show_col_types = FALSE)
shared_cols2 <- intersect(names(df2), names(df2_alff))
df2_combined <- bind_rows(df2 %>% select(all_of(shared_cols2)),
                          df2_alff %>% select(all_of(shared_cols2)))
cat("=== RESIDUALIZED MODEL: families present after merge ===\n")
print(table(df2_combined$model_family))

tbl2 <- df2_combined %>%
  filter(term %in% c("behav_resid", "behav_resid:groupstroke")) %>%
  filter(p.value < p_threshold) %>%
  mutate(family_order_idx = match(model_family, family_order)) %>%
  arrange(family_order_idx, feature_type, p.value) %>%
  mutate(Family  = family_label(model_family),
         Feature = clean_feature(feature_name),
         Term    = clean_term(term),
         Beta    = fmt_beta(estimate),
         SE      = fmt_se(std.error),
         t       = fmt_t(statistic),
         p       = fmt_p(p.value),
         `FDR q` = fmt_q(p_fdr),
         Sig     = sig_star(p_fdr, p.value),
         n_obs   = n_obs,
         n_subj  = n_subjects) %>%
  select(Family, Feature, Term, Beta, SE, t, p, `FDR q`, Sig, n_obs, n_subj)

df3_fc   <- read_csv(chaco_fc_path,   show_col_types = FALSE)
df3_alff <- read_csv(chaco_alff_path, show_col_types = FALSE)
shared_cols3 <- intersect(names(df3_fc), names(df3_alff))
df3_combined <- bind_rows(df3_fc   %>% select(all_of(shared_cols3)),
                          df3_alff %>% select(all_of(shared_cols3)))
cat("=== CHACO MODERATION: families present after merge ===\n")
print(table(df3_combined$model_family))

chaco_terms <- df3_combined %>%
  filter(term == "behav_z:chaco_z") %>%
  filter(model_family %in% c("FC_SEED", "FC_NETWORK")) %>%
  mutate(family_order_idx = match(model_family, family_order)) %>%
  arrange(family_order_idx, feature_type, p.value) %>%
  mutate(Family  = family_label(model_family),
         Feature = clean_feature(feature_name),
         Term    = "Vitality x ChaCo (interaction)",
         Beta    = fmt_beta(estimate),
         SE      = fmt_se(std.error),
         t       = fmt_t(statistic),
         p       = fmt_p(p.value),
         `FDR q` = fmt_q(p_fdr),
         Sig     = sig_star(p_fdr, p.value),
         n_obs   = n_obs,
         n_subj  = n_subjects) %>%
  select(Family, Feature, Term, Beta, SE, t, p, `FDR q`, Sig, n_obs, n_subj)

# =============================================================================
#  PROCESS TABLES S4 & S5 — Group differences (all metrics in one file)
# =============================================================================
df_group <- read_csv(alff_falff_path, show_col_types = FALSE)

cat("=== GROUP DIFF: model families present ===\n")
print(table(df_group$model_family))
cat(sprintf("sig_unc rows: %d\n", sum(df_group$sig_unc == TRUE)))

# S4 — FC Node Strength
tbl4 <- df_group %>%
  filter(model_family == "FC_NODE_STRENGTH") %>%
  filter(sig_unc == TRUE) %>%
  arrange(p_value) %>%
  mutate(
    Metric  = "FC Node Strength",
    Region  = clean_feature(parcel),
    t       = fmt_t(t_value),
    p       = fmt_p(p_value),
    `FDR q` = fmt_q(p_fdr),
    Sig     = sig_star(p_fdr, p_value),
    Est     = fmt_beta(estimate),
    SE      = fmt_se(std.error),
    n_obs   = n_obs,
    n_subj  = n_subjects
  ) %>%
  select(Metric, Region, t, p, `FDR q`, Sig, Est, SE, n_obs, n_subj)

# S5 — ALFF/fALFF
tbl5 <- df_group %>%
  filter(model_family %in% c("ALFF_ROI", "FALFF_ROI")) %>%
  filter(sig_unc == TRUE) %>%
  mutate(metric_idx = match(model_family, c("ALFF_ROI", "FALFF_ROI"))) %>%
  arrange(metric_idx, p_value) %>%
  mutate(
    Metric  = case_when(
      model_family == "ALFF_ROI"  ~ "ALFF",
      model_family == "FALFF_ROI" ~ "fALFF"
    ),
    Region  = clean_feature(parcel),
    t       = fmt_t(t_value),
    p       = fmt_p(p_value),
    `FDR q` = fmt_q(p_fdr),
    Sig     = sig_star(p_fdr, p_value),
    Est     = fmt_beta(estimate),
    SE      = fmt_se(std.error),
    n_obs   = n_obs,
    n_subj  = n_subjects
  ) %>%
  select(Metric, Region, t, p, `FDR q`, Sig, Est, SE, n_obs, n_subj)

# =============================================================================
#  BUILD WORKBOOK AND WRITE ALL SHEETS
# =============================================================================
wb  <- createWorkbook()
sty <- make_styles()

write_table_sheet(
  wb, "Table S1 - Primary model", tbl1,
  "Table S1. Vitality-imaging associations: depression-adjusted primary model (n = 201 observations, 132 subjects)",
  paste0("Models: Imaging Feature(ij) = b0 + b1*Vitality(ij) + b2*GDSS(ij) + b3*Group(i) + ",
         "b4*(Vitality x Group)(ij) + b5*Age(i) + b6*Sex(i) + b7*Time(ij) + u(i) + e(ij). ",
         "FDR correction applied within model family (Benjamini-Hochberg). ",
         "* FDR q < 0.05; + uncorrected p < 0.05 only. ",
         "Vitality scored such that lower values = greater fatigue. GDSS = Geriatric Depression Scale Score."),
  sty
)

write_table_sheet(
  wb, "Table S2 - Residualized model", tbl2,
  "Table S2. Vitality-imaging associations: residualized vitality sensitivity analysis (n = 201 observations, 132 subjects)",
  paste0("Vitality residualized with respect to GDSS and covariates: ",
         "Vitality(resid) = residuals(Vitality ~ GDSS + Age + Sex + Group + Time). ",
         "Models: Imaging Feature(ij) = b0 + b1*Vitality_resid(ij) + b2*GDSS(ij) + b3*Group(i) + ",
         "b4*(Vitality_resid x Group)(ij) + b5*Age(i) + b6*Sex(i) + b7*Time(ij) + u(i) + e(ij). ",
         "FDR correction applied within model family (Benjamini-Hochberg). ",
         "* FDR q < 0.05; + uncorrected p < 0.05 only."),
  sty
)

write_table_sheet(
  wb, "Table S3 - ChaCo moderation", chaco_terms,
  "Table S3. Structural disconnection moderation: Vitality x ChaCo interaction (stroke participants only)",
  paste0("Stroke-only moderation models restricted to between-network FC and ROI-based FC families. ",
         "Vitality, ChaCo, and GDSS z-scored within stroke sample. ",
         "Model: Imaging Feature(ij) = b0 + b1*behav_z(ij) + b2*chaco_z(i) + b3*gdss_z(ij) + ",
         "b4*(behav_z x chaco_z)(ij) + b5*Age(i) + b6*Sex(i) + b7*Time(ij) + u(i) + e(ij). ",
         "Only the behav_z x chaco_z interaction term is shown. ",
         "FDR correction applied within model family (Benjamini-Hochberg). ",
         "* FDR q < 0.05; + uncorrected p < 0.05 only."),
  sty
)

write_group_diff_sheet(
  wb, "Table S4 - FC Node Strength", tbl4,
  "Table S4. FC node strength group differences: stroke vs. control (GDSS-adjusted, uncorrected p < 0.05)",
  paste0("Group differences in regional ALFF and fALFF adjusted for depressive symptoms (GDSS), age, sex, and time. ",
         "Positive t-values indicate greater amplitude in stroke relative to controls. ",
         "No regions survived FDR correction. + uncorrected p < 0.05 only."),
  sty
)

write_group_diff_sheet(
  wb, "Table S5 - ALFF fALFF", tbl5,
  "Table S5. ALFF and fALFF group differences: stroke vs. control (GDSS-adjusted, uncorrected p < 0.05)",
  paste0("Group differences in parcel-wise FC node strength (sum of Fisher-z transformed correlations) ",
         "adjusted for depressive symptoms (GDSS), age, sex, and time. ",
         "Positive t-values indicate greater node strength in stroke relative to controls. ",
         "No regions survived FDR correction. + uncorrected p < 0.05 only."),
  sty
)

# =============================================================================
#  SAVE
# =============================================================================
saveWorkbook(wb, out_path, overwrite = TRUE)
cat(sprintf("\u2705 Saved: %s\n", out_path))

cat("\n=== SUMMARY ===\n")
for (nm in list(list("S1 Primary",     tbl1),
                list("S2 Residualized", tbl2),
                list("S3 ChaCo",        chaco_terms),
                list("S4 FC Node",      tbl4),
                list("S5 ALFF/fALFF",  tbl5))) {
  tbl <- nm[[2]]
  cat(sprintf("Table %s: %d rows | FDR sig: %d | Unc p<0.05: %d\n",
              nm[[1]], nrow(tbl),
              sum(tbl$Sig == "*"),
              sum(tbl$Sig == "+")))
}
