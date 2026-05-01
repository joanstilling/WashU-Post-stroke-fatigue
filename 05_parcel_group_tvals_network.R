suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lme4)
  library(broom.mixed)
  library(ggplot2)
  library(tidyr)
})

# ============================================================
# Paths
# ============================================================
merged_all_path <- "/Users/.../merged_imaging_behaviour_long_ALL.csv"
out_dir         <- "/Users/.../figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat_fc  <- read_csv(merged_fc_path, show_col_types = FALSE)
dat_all <- read_csv(merged_all_path, show_col_types = FALSE)

# ============================================================
# Combine features of interest
# ============================================================
dat <- bind_rows(
  dat_all %>% filter(feature_type %in% c("ALFF_ROI", "FALFF_ROI")),
  dat_fc  %>% filter(feature_type %in% c("FC_NODE_STRENGTH"))
)

# ============================================================
# Exact fs86 ROI -> Yeo9 mapping from your MATLAB code
# ============================================================
region_names <- c(
  "Left-Cerebellum-Cortex", "Left-Thalamus-Proper", "Left-Caudate",
  "Left-Putamen", "Left-Pallidum", "Left-Hippocampus", "Left-Amygdala",
  "Left-Accumbens-area", "Left-VentralDC", "Right-Cerebellum-Cortex",
  "Right-Thalamus-Proper", "Right-Caudate", "Right-Putamen",
  "Right-Pallidum", "Right-Hippocampus", "Right-Amygdala",
  "Right-Accumbens-area", "Right-VentralDC", "ctx-lh-bankssts",
  "ctx-lh-caudalanteriorcingulate", "ctx-lh-caudalmiddlefrontal",
  "ctx-lh-cuneus", "ctx-lh-entorhinal", "ctx-lh-fusiform",
  "ctx-lh-inferiorparietal", "ctx-lh-inferiortemporal",
  "ctx-lh-isthmuscingulate", "ctx-lh-lateraloccipital",
  "ctx-lh-lateralorbitofrontal", "ctx-lh-lingual",
  "ctx-lh-medialorbitofrontal", "ctx-lh-middletemporal",
  "ctx-lh-parahippocampal", "ctx-lh-paracentral",
  "ctx-lh-parsopercularis", "ctx-lh-parsorbitalis",
  "ctx-lh-parstriangularis", "ctx-lh-pericalcarine",
  "ctx-lh-postcentral", "ctx-lh-posteriorcingulate",
  "ctx-lh-precentral", "ctx-lh-precuneus",
  "ctx-lh-rostralanteriorcingulate", "ctx-lh-rostralmiddlefrontal",
  "ctx-lh-superiorfrontal", "ctx-lh-superiorparietal",
  "ctx-lh-superiortemporal", "ctx-lh-supramarginal", "ctx-lh-frontalpole",
  "ctx-lh-temporalpole", "ctx-lh-transversetemporal", "ctx-lh-insula",
  "ctx-rh-bankssts", "ctx-rh-caudalanteriorcingulate",
  "ctx-rh-caudalmiddlefrontal", "ctx-rh-cuneus", "ctx-rh-entorhinal",
  "ctx-rh-fusiform", "ctx-rh-inferiorparietal",
  "ctx-rh-inferiortemporal", "ctx-rh-isthmuscingulate",
  "ctx-rh-lateraloccipital", "ctx-rh-lateralorbitofrontal",
  "ctx-rh-lingual", "ctx-rh-medialorbitofrontal",
  "ctx-rh-middletemporal", "ctx-rh-parahippocampal", "ctx-rh-paracentral",
  "ctx-rh-parsopercularis", "ctx-rh-parsorbitalis",
  "ctx-rh-parstriangularis", "ctx-rh-pericalcarine", "ctx-rh-postcentral",
  "ctx-rh-posteriorcingulate", "ctx-rh-precentral", "ctx-rh-precuneus",
  "ctx-rh-rostralanteriorcingulate", "ctx-rh-rostralmiddlefrontal",
  "ctx-rh-superiorfrontal", "ctx-rh-superiorparietal",
  "ctx-rh-superiortemporal", "ctx-rh-supramarginal", "ctx-rh-frontalpole",
  "ctx-rh-temporalpole", "ctx-rh-transversetemporal", "ctx-rh-insula"
)

network_labels <- c(
  9,8,8,8,8,8,8,8,8,9,8,8,8,8,8,8,8,8,
  7,4,7,1,5,1,7,5,7,1,5,1,5,7,1,2,7,7,7,
  1,2,4,2,7,7,6,7,3,2,4,5,5,2,4,2,4,6,1,5,
  1,7,5,7,1,5,1,5,7,1,2,6,7,6,1,2,4,2,7,7,
  6,7,3,2,4,5,5,2,4
)

network_names <- c("VIS","SOM","DAN","VAN","LIM","FP","DMN","SUB","CER")

roi_lookup <- tibble(
  region_name = region_names,
  network_label = network_labels
) %>%
  mutate(
    network = network_names[network_label]
  )

# ============================================================
# Helpers for parcel name cleaning
# ============================================================
clean_parcel <- function(x) {
  x %>%
    str_remove("^NS_") %>%
    str_remove("^ALFF_") %>%
    str_remove("^fALFF_")
}

extract_hemi <- function(x) {
  x2 <- clean_parcel(x)
  case_when(
    str_detect(x2, "^ctx-lh-") ~ "LH",
    str_detect(x2, "^ctx-rh-") ~ "RH",
    str_detect(x2, "^Left-")   ~ "LH",
    str_detect(x2, "^Right-")  ~ "RH",
    TRUE ~ "Other"
  )
}

# ============================================================
# Prepare network-level data
# ============================================================
dat_net <- dat %>%
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
    group = factor(group, levels = c("control", "stroke")),
    sex   = factor(sex, levels = c("female", "male")),
    age   = suppressWarnings(as.numeric(age)),
    GDSS  = suppressWarnings(as.numeric(gdss_score)),
    time_mo = suppressWarnings(as.numeric(time_mo)),
    value = suppressWarnings(as.numeric(value)),
    study_id = as.character(study_id),
    event_name = as.character(event_name),
    parcel_clean = clean_parcel(feature_name),
    hemi = extract_hemi(feature_name)
  ) %>%
  left_join(
    roi_lookup %>% select(region_name, network),
    by = c("parcel_clean" = "region_name")
  ) %>%
  filter(!is.na(network)) %>%
  group_by(study_id, event_name, time_mo, group, age, sex, GDSS, feature_type, network) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    n_parcels = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    network = factor(network, levels = c("VIS","SOM","DAN","VAN","LIM","FP","DMN","SUB","CER")),
    feature_type = factor(feature_type, levels = c("ALFF_ROI","FALFF_ROI","FC_NODE_STRENGTH"))
  )

write_csv(dat_net, file.path(out_dir, "NETWORK_LEVEL_feature_means_with_GDSS.csv"))

# ============================================================
# Mixed model with GDSS covariate
# ============================================================
fit_network_model <- function(df_net_one) {
  
  df_net_one <- df_net_one %>%
    filter(!is.na(value), !is.na(group), !is.na(age), !is.na(sex), !is.na(study_id), !is.na(GDSS))
  
  if (nrow(df_net_one) < 12) return(NULL)
  if (n_distinct(df_net_one$group) < 2) return(NULL)
  if (n_distinct(df_net_one$study_id) < 6) return(NULL)
  
  use_time_mo <- sum(!is.na(df_net_one$time_mo)) > 0
  
  form <- if (use_time_mo) {
    value ~ group + age + sex + GDSS + time_mo + (1 | study_id)
  } else {
    value ~ group + age + sex + GDSS + factor(event_name) + (1 | study_id)
  }
  
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m <- try(lmer(form, data = df_net_one, REML = FALSE, control = ctrl), silent = TRUE)
  
  if (inherits(m, "try-error")) return(NULL)
  
  td <- broom.mixed::tidy(m, effects = "fixed")
  
  out <- td %>%
    filter(term == "groupstroke") %>%
    transmute(
      estimate   = estimate,
      std.error  = std.error,
      t_value    = statistic,
      p_value    = 2 * pnorm(abs(statistic), lower.tail = FALSE),
      n_obs      = nrow(df_net_one),
      n_subjects = n_distinct(df_net_one$study_id),
      mean_GDSS  = mean(df_net_one$GDSS, na.rm = TRUE)
    )
  
  if (nrow(out) == 0) return(NULL)
  out
}

# ============================================================
# Run feature x network models
# ============================================================
net_results <- dat_net %>%
  group_by(feature_type, network) %>%
  group_modify(~{
    out <- fit_network_model(.x)
    if (is.null(out)) return(tibble())
    out
  }) %>%
  ungroup() %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "fdr"),
    sig_unc = p_value < 0.05,
    sig_fdr = p_fdr < 0.05
  ) %>%
  arrange(feature_type, network)

write_csv(net_results, file.path(out_dir, "NETWORK_LEVEL_group_tvals_with_GDSS.csv"))
print(net_results)

# ============================================================
# Heatmap
# ============================================================
heat_df <- net_results %>%
  mutate(
    feature_type = factor(feature_type, levels = c("ALFF_ROI","FALFF_ROI","FC_NODE_STRENGTH")),
    network = factor(network, levels = c("VIS","SOM","DAN","VAN","LIM","FP","DMN","SUB","CER"))
  )

p_net <- ggplot(heat_df, aes(x = network, y = feature_type, fill = t_value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(
    data = heat_df %>% filter(sig_unc),
    fill = NA, color = "#00A651", linewidth = 1.0
  ) +
  geom_point(
    data = heat_df %>% filter(sig_fdr),
    shape = 8, size = 3, color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    midpoint = 0,
    name = "t-value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11)
  ) +
  labs(
    title = "Network-level stroke vs control effects across all timepoints",
    subtitle = "Adjusted for age, sex, GDSS, and time; green outline = p < 0.05; black star = FDR q < 0.05"
  )

print(p_net)

ggsave(
  filename = file.path(out_dir, "NETWORK_LEVEL_ALFF_FALFF_FC_heatmap_with_GDSS.png"),
  plot = p_net,
  width = 11,
  height = 4.8,
  dpi = 300
)

# ============================================================
# Wide summary table
# ============================================================
net_summary_wide <- net_results %>%
  select(feature_type, network, t_value, p_value, p_fdr) %>%
  mutate(
    summary = sprintf("t=%.2f, p=%.3f, q=%.3f", t_value, p_value, p_fdr)
  ) %>%
  select(feature_type, network, summary) %>%
  pivot_wider(names_from = network, values_from = summary)

write_csv(net_summary_wide, file.path(out_dir, "NETWORK_LEVEL_summary_wide_with_GDSS.csv"))
print(net_summary_wide)