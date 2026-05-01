suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lme4)
  library(broom.mixed)
})

merged_all_path <- "/Users/.../merged_imaging_behaviour_long_ALL.csv"
out_dir         <- "/Users/..."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat_fc  <- readr::read_csv(merged_fc_path,  show_col_types = FALSE)
dat_all <- readr::read_csv(merged_all_path, show_col_types = FALSE)

cat("dat_fc rows:",  nrow(dat_fc),  " | feature types:", paste(unique(dat_fc$feature_type),  collapse = ", "), "\n")
cat("dat_all rows:", nrow(dat_all), " | feature types:", paste(unique(dat_all$feature_type), collapse = ", "), "\n\n")

# ------------------------------------------------------------
export_parcel_group_tvals <- function(dat, feature_type_use, out_dir) {
  
  message("Running feature_type: ", feature_type_use)
  
  # Debug: check raw values before recode
  raw_group_vals <- sort(unique(as.character(dat$group)))
  message("  Raw group values in input data: ",
          paste(head(raw_group_vals, 10), collapse = ", "))
  
  df <- dat %>%
    mutate(
      group = case_when(
        group %in% c("0", 0, "stroke",  "Stroke")  ~ "stroke",
        group %in% c("1", 1, "control", "Control") ~ "control",
        TRUE ~ as.character(group)
      ),
      sex = case_when(
        sex %in% c("0", 0, "female", "Female") ~ "female",
        sex %in% c("1", 1, "male",   "Male")   ~ "male",
        TRUE ~ as.character(sex)
      ),
      group      = factor(group, levels = c("control", "stroke")),
      sex        = factor(sex,   levels = c("female",  "male")),
      age        = suppressWarnings(as.numeric(age)),
      value      = suppressWarnings(as.numeric(value)),
      time_mo    = suppressWarnings(as.numeric(time_mo)),
      gdss_score = suppressWarnings(as.numeric(gdss_score)),
      parcel     = as.character(feature_name),
      study_id   = as.character(study_id),
      event_name = as.character(event_name)
    )
  
  message("  Rows after mutate: ", nrow(df))
  message("  Rows with feature_type == '", feature_type_use, "': ",
          sum(df$feature_type == feature_type_use, na.rm = TRUE))
  
  df <- df %>%
    filter(feature_type == feature_type_use)
  
  message("  Rows after feature_type filter: ", nrow(df))
  
  df <- df %>%
    filter(!is.na(value), !is.na(group), !is.na(age),
           !is.na(sex), !is.na(study_id), !is.na(gdss_score))
  
  message("  Rows after complete-case filter: ", nrow(df))
  
  if (nrow(df) == 0) {
    stop("No rows found for feature_type: ", feature_type_use)
  }
  
  use_time_mo <- sum(!is.na(df$time_mo)) > 0
  
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  
  fit_one_parcel <- function(dfp) {
    if (n_distinct(dfp$group) < 2) return(NULL)
    if (nrow(dfp) < 12)            return(NULL)
    

    form <- if (use_time_mo) {
      value ~ group + age + sex + gdss_score + time_mo + (1 | study_id)
    } else {
      value ~ group + age + sex + gdss_score + factor(event_name) + (1 | study_id)
    }
    
    m <- try(lmer(form, data = dfp, REML = FALSE, control = ctrl),
             silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    
    td <- broom.mixed::tidy(m, effects = "fixed")
    
    out <- td %>%
      filter(term == "groupstroke") %>%
      transmute(
        t_value    = statistic,
        p_value    = 2 * pnorm(abs(statistic), lower.tail = FALSE),
        estimate   = estimate,
        std.error  = std.error,
        n_obs      = nrow(dfp),
        n_subjects = n_distinct(dfp$study_id)
      )
    
    if (nrow(out) == 0) return(NULL)
    out
  }
  
  parcel_stats <- df %>%
    group_by(parcel) %>%
    group_modify(~{
      out <- fit_one_parcel(.x)
      if (is.null(out)) return(tibble())
      out %>% mutate(parcel = unique(.x$parcel))
    }) %>%
    ungroup() %>%
    mutate(
      hemi = case_when(
        str_detect(parcel, "^ctx-lh-") ~ "lh",
        str_detect(parcel, "^ctx-rh-") ~ "rh",
        TRUE ~ "subcort"
      ),
     
      p_fdr    = p.adjust(p_value, method = "BH"),
      sig_unc  = p_value < 0.05,
      sig_fdr  = p_fdr   < 0.05,
      model_family = feature_type_use
    ) %>%
    arrange(p_value)
  
  # Diagnostic summary
  message(sprintf("  Total parcels: %d", nrow(parcel_stats)))
  message(sprintf("  Uncorrected p < 0.05: %d", sum(parcel_stats$sig_unc, na.rm = TRUE)))
  message(sprintf("  FDR q < 0.05: %d",          sum(parcel_stats$sig_fdr,  na.rm = TRUE)))
  

  out_csv <- file.path(out_dir,
                       paste0("ALLTIME_", feature_type_use,
                              "_group_tvals_GDSSadj.csv"))
  write_csv(
    parcel_stats %>%
      select(parcel, model_family, t_value, p_value, p_fdr,
             sig_unc, sig_fdr, hemi,
             estimate, std.error, n_obs, n_subjects),
    out_csv
  )
  
  message("Wrote: ", out_csv)
  print(head(parcel_stats, 10))
  
  invisible(parcel_stats)
}

# ------------------------------------------------------------
# Run for all three parcel-wise feature types
# Output files now include _GDSSadj suffix to distinguish
# from original non-adjusted versions
# ------------------------------------------------------------
cat("=== GDSS CHECK IN dat_all ===\n")
cat("Column 'gdss_score' exists:", "gdss_score" %in% names(dat_all), "\n")
if ("gdss_score" %in% names(dat_all)) {
  cat("Non-missing gdss_score rows:", sum(!is.na(dat_all$gdss_score)), "of", nrow(dat_all), "\n")
} else {
  cat("Available columns with 'gdss' or 'gds':",
      paste(grep("gdss|gds", names(dat_all), value = TRUE, ignore.case = TRUE),
            collapse = ", "), "\n")
}

# Check GDSS coverage specifically within ALFF_ROI rows
alff_check <- dat_all %>% filter(feature_type == "ALFF_ROI")
cat("\nALFF_ROI rows:", nrow(alff_check), "\n")
cat("ALFF_ROI rows with non-missing gdss_score:",
    sum(!is.na(alff_check$gdss_score)), "\n")
cat("ALFF_ROI rows with non-missing group, age, sex, value, gdss_score:",
    sum(!is.na(alff_check$gdss_score) & !is.na(alff_check$group) &
          !is.na(alff_check$age) & !is.na(alff_check$sex) &
          !is.na(alff_check$value)), "\n\n")

alff_stats  <- export_parcel_group_tvals(dat_all, "ALFF_ROI",         out_dir)
falff_stats <- export_parcel_group_tvals(dat_all, "FALFF_ROI",        out_dir)
fcns_stats  <- export_parcel_group_tvals(dat_fc,  "FC_NODE_STRENGTH",  out_dir)

# ------------------------------------------------------------
# Combined summary across all three families
# ------------------------------------------------------------
all_stats <- bind_rows(alff_stats, falff_stats, fcns_stats)

write_csv(
  all_stats %>%
    select(model_family, parcel, hemi, t_value, p_value, p_fdr,
           sig_unc, sig_fdr, estimate, std.error, n_obs, n_subjects),
  file.path(out_dir, "ALLTIME_ALL_group_tvals_GDSSadj_combined.csv")
)

cat("\n=== SUMMARY ACROSS ALL FAMILIES ===\n")
all_stats %>%
  group_by(model_family) %>%
  summarise(
    n_parcels   = n(),
    n_unc       = sum(sig_unc, na.rm = TRUE),
    n_fdr       = sum(sig_fdr, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  print()


# ============================================================
# Surface visualization 
# ============================================================
suppressPackageStartupMessages({
  library(ggseg)
  library(ggplot2)
  library(sf)
})

plot_surface_map <- function(feature_type_use, out_dir,
                             show_fdr_outline = TRUE) {
  
  in_csv <- file.path(out_dir,
                      paste0("ALLTIME_", feature_type_use,
                             "_group_tvals_GDSSadj.csv"))
  
  stats <- readr::read_csv(in_csv, show_col_types = FALSE)
  
  plot_df <- stats %>%
    mutate(
      t_value = as.numeric(t_value),
      p_value = as.numeric(p_value),
      p_fdr   = as.numeric(p_fdr)
    ) %>%
    filter(grepl("^ctx-[lr]h-|^NS_ctx-[lr]h-", parcel)) %>%
    mutate(
      hemi   = ifelse(grepl("lh", parcel), "left", "right"),
      region = parcel %>%
        str_remove("^NS_") %>%
        str_remove("^ctx-lh-") %>%
        str_remove("^ctx-rh-") %>%
        dplyr::recode(
          bankssts                 = "banks sts",
          caudalanteriorcingulate  = "caudal anterior cingulate",
          caudalmiddlefrontal      = "caudal middle frontal",
          frontalpole              = "frontal pole",
          inferiorparietal         = "inferior parietal",
          inferiortemporal         = "inferior temporal",
          isthmuscingulate         = "isthmus cingulate",
          lateraloccipital         = "lateral occipital",
          lateralorbitofrontal     = "lateral orbitofrontal",
          medialorbitofrontal      = "medial orbitofrontal",
          middletemporal           = "middle temporal",
          parahippocampal          = "parahippocampal",
          parsopercularis          = "pars opercularis",
          parsorbitalis            = "pars orbitalis",
          parstriangularis         = "pars triangularis",
          pericalcarine            = "pericalcarine",
          posteriorcingulate       = "posterior cingulate",
          rostralanteriorcingulate = "rostral anterior cingulate",
          rostralmiddlefrontal     = "rostral middle frontal",
          superiorfrontal          = "superior frontal",
          superiorparietal         = "superior parietal",
          superiortemporal         = "superior temporal",
          temporalpole             = "temporal pole",
          transversetemporal       = "transverse temporal",
          .default = .
        )
    )
  
  atlas_joined <- brain_join(plot_df, dk(), by = c("hemi", "region")) %>%
    filter(!is.na(parcel), !sf::st_is_empty(geometry))
  
  # Uncorrected outline
  sig_unc_joined <- atlas_joined %>%
    filter(!is.na(p_value) & p_value < 0.05, !sf::st_is_empty(geometry))
  
  # FDR outline
  sig_fdr_joined <- atlas_joined %>%
    filter(!is.na(p_fdr) & p_fdr < 0.05, !sf::st_is_empty(geometry))
  
  message(sprintf("%s: %d parcels mapped, %d unc sig, %d FDR sig",
                  feature_type_use,
                  nrow(atlas_joined),
                  nrow(sig_unc_joined),
                  nrow(sig_fdr_joined)))
  
  p <- ggplot() +
    geom_sf(data  = atlas_joined,
            aes(fill = t_value),
            color = "grey75", linewidth = 0.2) +
    scale_fill_gradient2(
      low      = "#3B4CC0",
      mid      = "white",
      high     = "#B40426",
      midpoint = 0,
      name     = "t-value",
      na.value = "grey90"
    ) +
    coord_sf(datum = NA, expand = FALSE) +
    theme_void() +
    labs(
      title    = sprintf("%s: stroke vs. control (GDSS-adjusted)",
                         feature_type_use),
      subtitle = "Green outline = uncorrected p < 0.05; black outline = FDR q < 0.05"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      legend.position = "right"
    )
  
  # Add uncorrected outline (green)
  if (nrow(sig_unc_joined) > 0) {
    p <- p + geom_sf(data      = sig_unc_joined,
                     fill      = NA,
                     color     = "chartreuse3",
                     linewidth = 0.9)
  }
  
  # Add FDR outline (black) if requested
  if (show_fdr_outline && nrow(sig_fdr_joined) > 0) {
    p <- p + geom_sf(data      = sig_fdr_joined,
                     fill      = NA,
                     color     = "black",
                     linewidth = 1.2)
  }
  
  out_png <- file.path(out_dir,
                       sprintf("%s_surface_map_GDSSadj.png",
                               feature_type_use))
  ggsave(out_png, p, width = 10, height = 6, dpi = 300, bg = "white")
  message("Saved: ", out_png)
  
  invisible(p)
}

# Run surface maps for all three feature types
plot_surface_map("ALFF_ROI",        out_dir)
plot_surface_map("FALFF_ROI",       out_dir)
plot_surface_map("FC_NODE_STRENGTH", out_dir)
