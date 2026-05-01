suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(lme4)
  library(cowplot)
})

# =============================================================================
#  PATHS  —  edit these 
# =============================================================================

paths <- list(
  fdr_main  = "/Users/.../ALL_terms_with_FDR_mixed.csv",
  fdr_resid = "/Users/.../ALL_vitality_resid_terms_with_FDR_mixed.csv",
  long_dat  = "/Users/.../merged_imaging_behaviour_long_ALL.csv",
  out_dir   = "/Users/..."
)

dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
#  SHARED HELPERS
# =============================================================================

recode_group_sex <- function(df) {
  df %>%
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
      group = factor(group, levels = c("control", "stroke")),
      sex   = factor(sex,   levels = c("female",  "male"))
    )
}

# Cast common numeric columns
cast_numeric <- function(df) {
  df %>%
    mutate(across(c(value, behav, gdss_score, age, time_mo), as.numeric))
}

# ggplot2 theme used across all panels
theme_vit <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size),
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
#  LOAD RESULTS TABLES
# =============================================================================

res_main  <- read_csv(paths$fdr_main,  show_col_types = FALSE)
res_resid <- read_csv(paths$fdr_resid, show_col_types = FALSE)

# =============================================================================
# FIGURE 1A — FP–LIM scatter with GDSS-adjusted regression lines
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
})

# -----------------------------------------------------------------------------
# OPTIONAL labels to print on the figure
# -----------------------------------------------------------------------------
p_fdr_label <- "FDR q = 0.048"   # edit to your final q
beta_label  <- "\u03b2 = 0.001" # edit to your final beta

# -----------------------------------------------------------------------------
# load and prepare long dataset
# -----------------------------------------------------------------------------
dat <- read_csv(paths$long_dat, show_col_types = FALSE) %>%
  recode_group_sex() %>%
  cast_numeric()

df_limpf <- dat %>%
  filter(feature_type == "NET_BETWEEN_FC", feature_name == "LIM-FP") %>%
  filter(
    !is.na(value), !is.na(behav), !is.na(gdss_score),
    !is.na(age), !is.na(sex), !is.na(time_mo), !is.na(study_id), !is.na(group)
  ) %>%
  droplevels() %>%
  mutate(
    group_plot = case_when(
      grepl("stroke",  tolower(as.character(group)))  ~ "Stroke",
      grepl("control", tolower(as.character(group)))  ~ "Control",
      TRUE ~ as.character(group)
    ),
    group_plot = factor(group_plot, levels = c("Control", "Stroke"))
  )

stopifnot(nrow(df_limpf) > 10)

# -----------------------------------------------------------------------------
# fit mixed model (matches analysis model)
# -----------------------------------------------------------------------------
m_limpf <- lmer(
  value ~ behav + group + behav * group + gdss_score + age + sex + time_mo + (1 | study_id),
  data = df_limpf,
  REML = FALSE
)

# Optional p-value extractor
get_coef_p <- function(model, term) {
  cc <- coef(summary(model))
  if (!term %in% rownames(cc)) return(NA_real_)
  if ("Pr(>|t|)" %in% colnames(cc)) return(cc[term, "Pr(>|t|)"])
  NA_real_
}

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}

p_vit <- get_coef_p(m_limpf, "behav")

# -----------------------------------------------------------------------------
# prediction grid: covariates fixed at median / mode
# -----------------------------------------------------------------------------
covariate_refs <- df_limpf %>%
  summarise(
    age        = median(age,        na.rm = TRUE),
    gdss_score = median(gdss_score, na.rm = TRUE),
    time_mo    = median(time_mo,    na.rm = TRUE),
    sex        = names(sort(table(sex), decreasing = TRUE))[1]
  )

grid_limpf <- crossing(
  covariate_refs,
  behav = seq(
    quantile(df_limpf$behav, 0.05, na.rm = TRUE),
    quantile(df_limpf$behav, 0.95, na.rm = TRUE),
    length.out = 100
  ),
  group = levels(df_limpf$group)
) %>%
  mutate(
    sex        = factor(sex, levels = levels(df_limpf$sex)),
    group      = factor(group, levels = levels(df_limpf$group)),
    group_plot = case_when(
      grepl("stroke",  tolower(as.character(group)))  ~ "Stroke",
      grepl("control", tolower(as.character(group)))  ~ "Control",
      TRUE ~ as.character(group)
    ),
    group_plot = factor(group_plot, levels = c("Control", "Stroke"))
  )

grid_limpf$pred <- predict(m_limpf, newdata = grid_limpf, re.form = NA)

# -----------------------------------------------------------------------------
# title/subtitle text
# -----------------------------------------------------------------------------
subtitle_main <- paste0(
  beta_label, "  •  ", p_fdr_label,
  "  •  n = ", nrow(df_limpf),
  "  •  subjects = ", dplyr::n_distinct(df_limpf$study_id),
  if (!is.na(p_vit)) paste0("  •  p = ", fmt_p(p_vit)) else ""
)

# -----------------------------------------------------------------------------
# shared style matching the other scatter figures
# -----------------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "black", hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(8, 8, 8, 8)
  )

shape_scale <- scale_shape_manual(
  values = c("Control" = 16, "Stroke" = 17),
  drop = FALSE
)

color_scale <- scale_color_manual(
  values = c("Control" = "grey45", "Stroke" = "black"),
  drop = FALSE
)

# -----------------------------------------------------------------------------
# plot
# -----------------------------------------------------------------------------
p_limpf <- ggplot(df_limpf, aes(x = behav, y = value)) +
  
  geom_point(
    aes(color = group_plot, shape = group_plot),
    alpha = 0.72,
    size = 2.2
  ) +
  
  geom_line(
    data = grid_limpf,
    aes(
      x = behav,
      y = pred,
      group = group_plot,
      linetype = group_plot
    ),
    color = "black",
    linewidth = 1.15
  ) +
  
  color_scale +
  shape_scale +
  
  scale_linetype_manual(
    values = c(
      "Control" = "solid",
      "Stroke"  = "dotted"
    )
  ) +
  
  labs(
    title = "Frontoparietal–limbic connectivity (depression adjusted)",
    subtitle = subtitle_main,
    x = "SF-36 Vitality (lower = greater fatigue)",
    y = "Connectivity (Fisher z)"
  ) +
  
  theme_pub

print(p_limpf)

# -----------------------------------------------------------------------------
# save
# -----------------------------------------------------------------------------
ggsave(
  file.path(paths$out_dir, "Figure1_LIM_FP_adjusted.tiff"),
  p_limpf,
  width = 5.8,
  height = 4.5,
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

ggsave(
  file.path(paths$out_dir, "Figure1_LIM_FP_adjusted.png"),
  p_limpf,
  width = 5.8,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(paths$out_dir, "Figure1_LIM_FP_adjusted.pdf"),
  p_limpf,
  width = 5.8,
  height = 4.5,
  bg = "white"
)

# =============================================================================
# FIGURE 1B — Between-network vitality matrix
# Triangular, grey diagonal, black = FDR, green = p<0.05, labeled key cells
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Canonical network order
net_order <- c("VIS", "SOM", "DAN", "VAN", "LIM", "FP", "DMN", "SUB", "CER")

# -------------------------------------------------------------------------
# Parse results
# -------------------------------------------------------------------------
net_mat_full <- res_main %>%
  filter(
    feature_type == "NET_BETWEEN_FC",
    effect_type  == "vitality_main"
  ) %>%
  separate(feature_name, into = c("net1", "net2"), sep = "-", remove = FALSE) %>%
  mutate(
    sig_fdr = !is.na(p_fdr) & p_fdr < 0.05,
    sig_p   = !is.na(p.value) & p.value < 0.05,
    net1 = factor(net1, levels = net_order),
    net2 = factor(net2, levels = net_order)
  )

# Add mirrored half
net_mat_mirror <- net_mat_full %>%
  transmute(
    feature_name,
    net1 = net2,
    net2 = net1,
    estimate,
    p.value,
    p_fdr,
    sig_fdr,
    sig_p
  )

net_mat_all <- bind_rows(
  net_mat_full %>%
    select(feature_name, net1, net2, estimate, p.value, p_fdr, sig_fdr, sig_p),
  net_mat_mirror
) %>%
  distinct(net1, net2, .keep_all = TRUE)

# -------------------------------------------------------------------------
# Full grid
# -------------------------------------------------------------------------
grid_df <- expand.grid(
  net1 = factor(net_order, levels = net_order),
  net2 = factor(net_order, levels = net_order)
) %>%
  as_tibble()

net_mat_all <- grid_df %>%
  left_join(net_mat_all, by = c("net1", "net2")) %>%
  mutate(
    i = as.integer(net1),
    j = as.integer(net2),
    is_diag = i == j,
    upper_triangle = j >= i
  )

# Keep upper triangle only
net_tri <- net_mat_all %>%
  filter(upper_triangle)

# Counts
n_sig_fdr <- net_mat_full %>% filter(sig_fdr) %>% nrow()
n_sig_p   <- net_mat_full %>% filter(sig_p & !sig_fdr) %>% nrow()

subtitle_txt <- paste0(
  "Black outlines: FDR q < 0.05 (n = ", n_sig_fdr,
  "); Green outlines: uncorrected p < 0.05 only (n = ", n_sig_p, ")"
)

# -------------------------------------------------------------------------
# LABEL THESE CELLS
# Put the pairs you want to call out here
# IMPORTANT: use canonical order as they appear in feature_name
# -------------------------------------------------------------------------
label_pairs <- tibble::tribble(
  ~net1, ~net2, ~label,
  "LIM", "FP",  "FP–LIM",
  "SUB", "VAN", "SUB–VAN"
  # add more only if truly necessary:
  # "SUB", "DMN", "SUB–DMN"
)

label_pairs <- label_pairs %>%
  mutate(
    net1 = factor(net1, levels = net_order),
    net2 = factor(net2, levels = net_order)
  )

# Because the displayed triangle uses reversed y-axis,
# we match labels to whichever orientation exists in net_tri
label_df <- net_tri %>%
  mutate(
    net1_chr = as.character(net1),
    net2_chr = as.character(net2)
  ) %>%
  inner_join(
    label_pairs %>%
      mutate(
        net1_chr = as.character(net1),
        net2_chr = as.character(net2)
      ) %>%
      select(net1_chr, net2_chr, label),
    by = c("net1_chr", "net2_chr")
  )

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
p_netmat_tri <- ggplot(net_tri, aes(x = net1, y = net2)) +
  
  # grey diagonal
  geom_tile(
    data = filter(net_tri, is_diag),
    fill = "grey85",
    colour = "white",
    linewidth = 0.35
  ) +
  
  # off-diagonal effect tiles
  geom_tile(
    data = filter(net_tri, !is_diag),
    aes(fill = estimate),
    colour = "white",
    linewidth = 0.35
  ) +
  
  # green outlines: nominal only
  geom_tile(
    data = filter(net_tri, !is_diag, sig_p & !sig_fdr),
    fill = NA,
    colour = "#1B9E77",
    linewidth = 0.9
  ) +
  
  # black outlines: FDR significant
  geom_tile(
    data = filter(net_tri, !is_diag, sig_fdr),
    fill = NA,
    colour = "black",
    linewidth = 1.05
  ) +
  
  # labels for key cells
  geom_text(
    data = label_df,
    aes(label = label),
    size = 2.5,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low      = "#3B4CC0",
    mid      = "white",
    high     = "#B40426",
    midpoint = 0,
    name     = expression(beta~"(Vitality [lower = greater fatigue])"),
    na.value = "grey95"
  ) +
  
  coord_fixed() +
  scale_y_discrete(limits = rev(net_order)) +
  
  labs(
    title = "Between-network connectivity associated with fatigue",
    subtitle = subtitle_txt,
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(p_netmat_tri)

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------
ggsave(
  file.path(paths$out_dir, "Fig_network_matrix_vitality_NET_BETWEEN_FC_triangular_labeled.png"),
  p_netmat_tri,
  width = 7.4,
  height = 6.1,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(paths$out_dir, "Fig_network_matrix_vitality_NET_BETWEEN_FC_triangular_labeled.tiff"),
  p_netmat_tri,
  width = 7.4,
  height = 6.1,
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

ggsave(
  file.path(paths$out_dir, "Fig_network_matrix_vitality_NET_BETWEEN_FC_triangular_labeled.pdf"),
  p_netmat_tri,
  width = 7.4,
  height = 6.1,
  bg = "white"
)

# =============================================================================
#  FIGURE 2 — Thalamic panels (residualised vitality)
# =============================================================================

# Residualise vitality against depression and covariates
# Use the full long dataset (one row per observation, ignoring imaging feature)
df_for_resid <- dat %>%
  distinct(study_id, time_mo, .keep_all = TRUE) %>%   # one row per person-timepoint
  filter(!is.na(behav), !is.na(gdss_score), !is.na(age),
         !is.na(sex),   !is.na(group),      !is.na(time_mo))

m_resid <- lm(
  behav ~ gdss_score + age + sex + group + time_mo,
  data = df_for_resid
)
df_for_resid$behav_resid <- resid(m_resid)

# Join residuals back to imaging data
df_imaging_resid <- dat %>%
  left_join(
    df_for_resid %>% select(study_id, time_mo, behav_resid),
    by = c("study_id", "time_mo")
  ) %>%
  filter(!is.na(behav_resid), !is.na(value))

# ---- Thalamus–Insula scatter ----
df_thal_ins <- df_imaging_resid %>%
  filter(feature_type == "SEED_EDGE_FC", feature_name == "THalamus<->Insula")

p_thal_ins <- ggplot(df_thal_ins, aes(x = behav_resid, y = value)) +
  geom_point(aes(shape = group), alpha = 0.65, size = 1.6) +
  geom_smooth(
    aes(linetype = group),
    method    = "lm",
    se        = FALSE,
    linewidth = 0.85
  ) +
  scale_shape_manual(values = c(control = 1, stroke = 16)) +
  labs(
    x        = "Vitality residual (VT | GDSS + covariates)",
    y        = "Thalamus–Insula connectivity (Fisher z)",
    shape    = NULL,
    linetype = NULL,
    title    = "Thalamus–Insula  (FDR q = 0.009)"
  ) +
  theme_vit()

# ---- Thalamus→VAN scatter ----
df_thal_van <- df_imaging_resid %>%
  filter(feature_type == "SEED_TO_NET_FC", feature_name == "SEED:THalamus->VAN")

p_thal_van <- ggplot(df_thal_van, aes(x = behav_resid, y = value)) +
  geom_point(aes(shape = group), alpha = 0.65, size = 1.6) +
  geom_smooth(
    aes(linetype = group),
    method    = "lm",
    se        = FALSE,
    linewidth = 0.85
  ) +
  scale_shape_manual(values = c(control = 1, stroke = 16)) +
  labs(
    x        = "Vitality residual (VT | GDSS + covariates)",
    y        = "Thalamus→VAN connectivity (Fisher z)",
    shape    = NULL,
    linetype = NULL,
    title    = "Thalamus→VAN  (FDR q = 0.010)"
  ) +
  theme_vit()

# =============================================================================
#  Assemble slope plot + thalamic panels
# =============================================================================

fig2_bottom <- plot_grid(
  p_thal_ins, p_thal_van,
  ncol   = 2,
  labels = c("B", "C"),
  label_size = 12
)

fig2 <- plot_grid(
  p_slope,
  fig2_bottom,
  ncol        = 1,
  rel_heights = c(1.2, 1),
  labels      = c("A", ""),
  label_size  = 12
)

ggsave(
  file.path(paths$out_dir, "Figure2_confirmatory_and_thalamic.tiff"),
  fig2,
  width = 7, height = 9, dpi = 600, compression = "lzw"
)

# =============================================================================
# FIGURE 3/S1: Surface maps — ALFF_ROI, FALFF_ROI, FC_NODE_STRENGTH
#  Source: ALLTIME_ALL_group_tvals_GDSSadj_combined.csv
#
#  2x2 grid per measure:
#    Row 1: Left lateral  | Right lateral
#    Row 2: Left medial   | Right medial
#  + Subcortical axial slices column
#  + Legend
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(sf); library(ggseg)
  library(patchwork); library(cowplot)
})

out_dir      <- "/Users/..."
# Single combined GDSS-adjusted file — all three feature types
combined_csv <- file.path(out_dir, "ALLTIME_ALL_group_tvals_GDSSadj_combined.csv")

# =============================================================================
#  1. LOAD ATLASES
# =============================================================================
dk_raw <- tryCatch(ggseg::dk, error = function(e) NULL)
if (is.function(dk_raw)) dk_raw <- ggseg::dk()
atlas_sf <- if (inherits(dk_raw, "brain_atlas")) dk_raw$data else as.data.frame(dk_raw)
atlas_sf  <- sf::st_as_sf(atlas_sf)
atlas_sf  <- atlas_sf[!sf::st_is_empty(atlas_sf), ]
message("Atlas view values: ", paste(sort(unique(atlas_sf$view)), collapse = ", "))

aseg_raw <- tryCatch(ggseg::aseg, error = function(e) NULL)
if (is.function(aseg_raw)) aseg_raw <- ggseg::aseg()
aseg_sf <- if (inherits(aseg_raw, "brain_atlas")) aseg_raw$data else as.data.frame(aseg_raw)
aseg_sf  <- sf::st_as_sf(aseg_sf)
aseg_sf  <- aseg_sf[!sf::st_is_empty(aseg_sf), ]
message("aseg regions: ", paste(sort(unique(aseg_sf$region)), collapse = ", "))

# Bounding boxes
get_bbox <- function(hemi_val, view_val) {
  sub <- atlas_sf[
    !is.na(atlas_sf$hemi) & atlas_sf$hemi == hemi_val &
      !is.na(atlas_sf$view) & atlas_sf$view == view_val, ]
  if (nrow(sub) == 0) return(NULL)
  bb    <- sf::st_bbox(sub)
  xpad  <- (bb["xmax"] - bb["xmin"]) * 0.03
  ypad  <- (bb["ymax"] - bb["ymin"]) * 0.03
  list(xlim = c(bb["xmin"] - xpad, bb["xmax"] + xpad),
       ylim = c(bb["ymin"] - ypad, bb["ymax"] + ypad))
}
bb <- list(
  left_lateral  = get_bbox("left",  "lateral"),
  right_lateral = get_bbox("right", "lateral"),
  left_medial   = get_bbox("left",  "medial"),
  right_medial  = get_bbox("right", "medial")
)

# =============================================================================
#  2. LOAD COMBINED STATS ONCE
# =============================================================================
combined_stats <- readr::read_csv(combined_csv, show_col_types = FALSE) %>%
  mutate(
    t_value = suppressWarnings(as.numeric(t_value)),
    p_value = suppressWarnings(as.numeric(p_value)),
    p_fdr   = suppressWarnings(as.numeric(p_fdr))
  )

message("Combined stats: ", nrow(combined_stats), " rows | families: ",
        paste(sort(unique(combined_stats$model_family)), collapse = ", "))

# =============================================================================
#  3. PARCEL NAME CLEANERS
# =============================================================================
recode_region <- function(x) {
  dplyr::recode(x,
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
                .default = x
  )
}

parse_parcels <- function(parcel_vec) {
  tibble(parcel = as.character(parcel_vec)) %>%
    mutate(
      pc     = parcel %>% str_remove("^NS_|^ALFF_|^fALFF_"),
      hemi   = case_when(
        str_detect(pc, "^ctx-lh-") ~ "left",
        str_detect(pc, "^ctx-rh-") ~ "right",
        TRUE ~ NA_character_
      ),
      region = recode_region(
        pc %>% str_remove("^ctx-lh-") %>% str_remove("^ctx-rh-"))
    ) %>%
    select(parcel, hemi, region)
}

parse_subcort <- function(parcel_vec) {
  tibble(parcel = as.character(parcel_vec)) %>%
    mutate(
      pc = parcel %>%
        str_remove("^NS_") %>% str_remove("^ALFF_") %>%
        str_remove("^fALFF_") %>% str_remove("^FC_"),
      region = pc %>%
        str_remove("^Left-") %>% str_remove("^Right-") %>%
        dplyr::recode(
          "Thalamus-Proper"   = "Thalamus Proper",
          "Caudate"           = "Caudate",
          "Putamen"           = "Putamen",
          "Pallidum"          = "Pallidum",
          "Hippocampus"       = "Hippocampus",
          "Amygdala"          = "Amygdala",
          "Accumbens-area"    = "accumbens area",
          "VentralDC"         = "ventraldc",
          "Cerebellum-Cortex" = "Cerebellum",
          .default = pc %>% str_remove("^Left-") %>% str_remove("^Right-")
        )
    ) %>%
    select(parcel, region)
}

# =============================================================================
#  4. JOIN FUNCTIONS — filter combined_stats by feature type
# =============================================================================
join_stats <- function(feature_type_use) {
  stats <- combined_stats %>% filter(model_family == feature_type_use)
  message("  ", feature_type_use, ": ", nrow(stats), " parcels")
  pm  <- parse_parcels(stats$parcel)
  df  <- stats[, !names(stats) %in% c("hemi", "region"), drop = FALSE]
  df$hemi   <- pm$hemi
  df$region <- pm$region
  df  <- df[!is.na(df$hemi) & !is.na(df$region), ]
  atlas_sf %>%
    left_join(df[, c("hemi", "region", "t_value", "p_value", "p_fdr")],
              by = c("hemi", "region")) %>%
    { .[!sf::st_is_empty(.), ] }
}

join_subcort <- function(feature_type_use) {
  stats <- combined_stats %>%
    filter(model_family == feature_type_use,
           str_detect(parcel, "Left-|Right-"))
  if (nrow(stats) == 0) return(NULL)
  pm  <- parse_subcort(stats$parcel)
  df  <- stats[, !names(stats) %in% "region", drop = FALSE]
  df$region <- pm$region
  df_avg <- df %>%
    group_by(region) %>%
    summarise(
      t_value = mean(t_value, na.rm = TRUE),
      p_value = mean(p_value, na.rm = TRUE),
      p_fdr   = mean(p_fdr,   na.rm = TRUE),
      .groups = "drop"
    )
  joined <- aseg_sf %>%
    left_join(df_avg[, c("region", "t_value", "p_value", "p_fdr")],
              by = "region")
  joined[!sf::st_is_empty(joined), ]
}

# =============================================================================
#  5. PANEL BUILDERS
# =============================================================================
make_panel <- function(ad, hemi_val, view_val, lims, sig_p = 0.05,
                       use_fdr = FALSE) {
  pd  <- ad[!is.na(ad$hemi) & ad$hemi == hemi_val &
              !is.na(ad$view) & ad$view == view_val, ]
  pd  <- pd[!sf::st_is_empty(pd), ]
  if (nrow(pd) == 0) return(ggplot() + theme_void())
  
  p_col <- if (use_fdr) "p_fdr" else "p_value"
  sig   <- pd[!is.na(pd[[p_col]]) & pd[[p_col]] < sig_p, ]
  
  bx_med   <- bb[[paste0(hemi_val, "_medial")]]
  bx_own   <- bb[[paste0(hemi_val, "_", view_val)]]
  med_w    <- diff(bx_med$xlim); med_h <- diff(bx_med$ylim)
  own_xmid <- mean(bx_own$xlim); own_ymid <- mean(bx_own$ylim)
  bx <- list(
    xlim = c(own_xmid - med_w/2, own_xmid + med_w/2),
    ylim = c(own_ymid - med_h/2, own_ymid + med_h/2)
  )
  
  p <- ggplot(pd) +
    geom_sf(aes(fill = t_value), colour = "grey78", linewidth = 0.15) +
    scale_fill_gradient2(
      low = "#2c7bb6", mid = "white", high = "#d7191c",
      midpoint = 0, limits = lims, oob = scales::squish,
      na.value = "grey92", name = "t-value", guide = "none"
    ) +
    coord_sf(datum = NA, expand = FALSE, xlim = bx$xlim, ylim = bx$ylim) +
    { if (hemi_val == "left" && view_val == "lateral")
      scale_x_reverse() else scale_x_continuous() } +
    theme_void() +
    theme(plot.margin = margin(3, 3, 3, 3))
  
  if (nrow(sig) > 0)
    p <- p + geom_sf(data = sig, fill = NA,
                     colour = "#00A651", linewidth = 2.5)
  p
}

make_subcort_panel <- function(ad_sub, lims, sig_p = 0.05,
                               use_fdr = FALSE) {
  if (is.null(ad_sub) || nrow(ad_sub) == 0) return(ggplot() + theme_void())
  
  p_col <- if (use_fdr) "p_fdr" else "p_value"
  sig   <- ad_sub[!is.na(ad_sub[[p_col]]) & ad_sub[[p_col]] < sig_p, ]
  ad_ax <- ad_sub[!is.na(ad_sub$view) & grepl("^axial", ad_sub$view), ]
  sig   <- sig[!is.na(sig$view) & grepl("^axial", sig$view), ]
  if (nrow(ad_ax) == 0) return(ggplot() + theme_void())
  
  p <- ggplot(ad_ax) +
    geom_sf(aes(fill = t_value), colour = "grey78", linewidth = 0.15) +
    scale_fill_gradient2(
      low = "#2c7bb6", mid = "white", high = "#d7191c",
      midpoint = 0, limits = lims, oob = scales::squish,
      na.value = "grey92", guide = "none"
    ) +
    facet_wrap(~view, nrow = 1) +
    coord_sf(datum = NA, expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(3,3,3,3), strip.text = element_blank())
  
  if (nrow(sig) > 0)
    p <- p + geom_sf(data = sig, fill = NA,
                     colour = "#00A651", linewidth = 1.2)
  p
}

# =============================================================================
#  6. BUILD GRID
# =============================================================================
make_legend <- function(lims) {
  df_leg <- data.frame(x=1, y=seq(lims[1],lims[2],length.out=100),
                       z=seq(lims[1],lims[2],length.out=100))
  p_leg <- ggplot(df_leg, aes(x=x, y=y, fill=z)) +
    geom_tile() +
    scale_fill_gradient2(
      low="#2c7bb6", mid="white", high="#d7191c",
      midpoint=0, limits=lims, name="t-value",
      guide=guide_colorbar(
        barheight=unit(5,"cm"), barwidth=unit(0.5,"cm"),
        title.position="top",
        title.theme=element_text(size=24, face="bold", fontfamily="arial"),
        label.theme=element_text(size=20)
      )
    ) +
    theme_void() + theme(legend.position="right")
  cowplot::get_legend(p_leg)
}

build_grid <- function(feature_type_use, label, lims, use_fdr = FALSE) {
  message("\n--- ", label, " ---")
  ad     <- join_stats(feature_type_use)
  ad_sub <- join_subcort(feature_type_use)
  
  p_ll <- make_panel(ad, "left",  "lateral", lims, use_fdr = use_fdr)
  p_rl <- make_panel(ad, "right", "lateral", lims, use_fdr = use_fdr)
  p_lm <- make_panel(ad, "left",  "medial",  lims, use_fdr = use_fdr)
  p_rm <- make_panel(ad, "right", "medial",  lims, use_fdr = use_fdr)
  
  leg <- make_legend(lims)
  
  title_p <- ggdraw() +
    draw_label(measure_labels[[label]], x=0.46, y=0.5,
               fontface="bold", size=36, fontfamily="arial", hjust=0.5)
  
  col_left  <- plot_grid(p_rl, p_rm, ncol=1, rel_heights=c(1,1))
  col_right <- plot_grid(p_ll, p_lm, ncol=1, rel_heights=c(1,1))
  
  # Subcortical: split axial slices into two rows to match cortical layout
  if (!is.null(ad_sub) && nrow(ad_sub) > 0) {
    axial_views <- sort(unique(ad_sub$view[grepl("^axial", ad_sub$view)]))
    split_at    <- ceiling(length(axial_views) / 2)
    views_top   <- axial_views[seq_len(split_at)]
    views_bot   <- axial_views[seq(split_at + 1, length(axial_views))]
    
    make_axial_stack <- function(views) {
      plots <- lapply(views, function(v) {
        slice <- ad_sub[!is.na(ad_sub$view) & ad_sub$view == v, ]
        slice <- slice[!sf::st_is_empty(slice), ]
        p_col <- if (use_fdr) "p_fdr" else "p_value"
        sig_s <- slice[!is.na(slice[[p_col]]) & slice[[p_col]] < 0.05, ]
        p <- ggplot(slice) +
          geom_sf(aes(fill = t_value), colour="grey78", linewidth=0.12) +
          scale_fill_gradient2(
            low="#2c7bb6", mid="white", high="#d7191c",
            midpoint=0, limits=lims, oob=scales::squish,
            na.value="grey92", guide="none") +
          coord_sf(datum=NA, expand=FALSE) +
          theme_void() + theme(plot.margin=margin(2,2,2,2))
        if (nrow(sig_s) > 0)
          p <- p + geom_sf(data=sig_s, fill=NA,
                           colour="#00A651", linewidth=1.2)
        p
      })
      plot_grid(plotlist=plots, ncol=1)
    }
    
    sub_col <- plot_grid(
      make_axial_stack(views_top),
      make_axial_stack(views_bot),
      ncol=1, rel_heights=c(1,1)
    )
  } else {
    sub_col <- ggplot() + theme_void()
  }
  
  brain_grid <- plot_grid(col_left, col_right, sub_col, leg,
                          nrow=1, rel_widths=c(2,2,0.7,0.4))
  plot_grid(title_p, brain_grid, ncol=1, rel_heights=c(0.06,1))
}

# =============================================================================
#  7. RUN
# =============================================================================
t_limits <- list(
  ALFF_ROI         = c(-2, 2),
  FALFF_ROI        = c(-3, 3),
  FC_NODE_STRENGTH = c(-2, 2)
)
measures <- c("ALFF_ROI", "FALFF_ROI", "FC_NODE_STRENGTH")

measure_labels <- c(
  ALFF_ROI         = "ALFF",
  FALFF_ROI        = "fALFF",
  FC_NODE_STRENGTH = "FC Node Strength"
)

# use_fdr = FALSE -> outline at uncorrected p < 0.05
# use_fdr = TRUE  -> outline at FDR q < 0.05
grids <- lapply(measures, function(ft) {
  tryCatch(
    build_grid(ft, ft, t_limits[[ft]], use_fdr = FALSE),
    error = function(e) { message("FAILED ", ft, ": ", e$message); NULL }
  )
})
names(grids) <- measures
grids <- Filter(Negate(is.null), grids)
if (length(grids) == 0) stop("No grids built.")

# =============================================================================
#  8. SAVE
# =============================================================================
for (nm in names(grids)) {
  ggsave(file.path(out_dir, paste0("SURFACE_grid_", nm, "_GDSSadj.png")),
         grids[[nm]], width=12, height=13, dpi=300, bg="white")
  ggsave(file.path(out_dir, paste0("SURFACE_grid_", nm, "_GDSSadj.tiff")),
         grids[[nm]], width=12, height=13, dpi=300,
         compression="lzw", bg="white")
  message("Saved: SURFACE_grid_", nm, "_GDSSadj")
}

combined <- plot_grid(
  plotlist         = grids,
  ncol             = 1,
  labels           = c("A", "B", "C"),
  label_size       = 28,
  label_fontface   = "bold",
  label_fontfamily = "Arial"
)
ggsave(file.path(out_dir, "SURFACE_grid_ALL_combined_GDSSadj.png"),
       combined, width=18, height=13*length(grids), dpi=300, bg="white")
ggsave(file.path(out_dir, "SURFACE_grid_ALL_combined_GDSSadj.tiff"),
       combined, width=18, height=13*length(grids), dpi=300,
       compression="lzw", bg="white")
message("Saved: combined. Done. Files in: ", out_dir)