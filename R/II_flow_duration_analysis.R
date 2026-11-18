#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Flow duration Analysis
# Date: 10/27/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Explore where samples were collected in flow duration space
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Approach
#  Step 1: Download all flow data
#  Step 2: Create flow duration curves
#  Step 3: Determine flow duration where samples were taken
#  Step 4: Create viz showing where samples curved flow duration curves

# 1.0 Setup workspace ----------------------------------------------------------
#Clear workspace
remove(list=ls())

#install libraries of interest
library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)

#Define data directory
data_dir <- "data/discharge"

# 2.0 Ingest and standardize flow data -----------------------------------------
# 2.1 Select files, preferring the wFlags version when one exists
q_files <- tibble(rel_path = list.files(data_dir, recursive = TRUE, pattern = "\\.csv$")) %>%
  filter(str_detect(rel_path, "_Q/")) %>%
  mutate(
    path     = file.path(data_dir, rel_path),
    region   = str_extract(rel_path, "^[A-Z]+"),
    file     = basename(rel_path),
    siteId   = str_extract(file, "[A-Z]{3}\\d{2}"),
    has_flag = str_detect(file, "wFlags"),
    file_key = str_remove(file, "_wFlags")
  ) %>%
  group_by(region, file_key) %>%
  slice_max(has_flag, n = 1, with_ties = FALSE) %>%
  ungroup()

q_files %>% select(region, file, has_flag) %>% print(n = Inf)

# 2.2 Read, standardize columns, and bind
read_q <- function(path, region, siteId, file) {
  df <- read_csv(path, col_types = cols(.default = "c"))
  
  # Add any columns that this file lacks
  for (col in c("datetime_UTC", "datetime", "wetdry", "QAQC", "na_flag")) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  
  df %>%
    mutate(
      region    = region,
      siteId    = siteId,
      file      = file,
      time_text = coalesce(datetime_UTC, datetime),
      # Pick the parser from the text pattern: ISO (2022-03-27...) or m/d/Y
      datetime  = if_else(
        str_detect(time_text, "^\\d{4}-\\d{2}-\\d{2}"),
        ymd_hms(time_text, tz = "UTC", truncated = 2, quiet = TRUE),
        mdy_hms(time_text, tz = "UTC", truncated = 2, quiet = TRUE)
      ),
      Q_Ls      = as.numeric(Q_Ls),
      # GBJ and JDR have no wetdry column, so build it from na_flag
      wetdry    = case_when(
        !is.na(wetdry)     ~ wetdry,
        na_flag %in% "Dry" ~ "dry",
        na_flag %in% "ND"  ~ NA_character_,
        !is.na(Q_Ls)       ~ "wet",
        TRUE               ~ NA_character_
      ),
      # Dry rows with missing discharge become zero flow (confirm in ReadMe)
      Q_Ls      = if_else(wetdry %in% "dry" & is.na(Q_Ls), 0, Q_Ls)
    ) %>%
    select(region, siteId, file, datetime, Q_Ls, wetdry, QAQC, na_flag)
}

q_sub <- q_files %>%
  select(path, region, siteId, file) %>%
  pmap_dfr(read_q)

# 2.3 Checks
# Rows that failed to parse or fall outside 2021 to 2024
q_sub %>%
  filter(is.na(datetime) | year(datetime) < 2021 | year(datetime) > 2024) %>%
  count(region, file, name = "n_bad_dates")

# Duplicate time stamps
q_sub %>%
  count(region, siteId, datetime) %>%
  filter(n > 1) %>%
  count(region, siteId, name = "n_duplicate_timestamps")

# Flag counts
q_sub %>% count(region, QAQC)   %>% print(n = Inf)
q_sub %>% count(region, na_flag) %>% print(n = Inf)

# 3.0 Aggregate to daily values, keeping flags ---------------------------------
# Flag severity, highest first. Unrecognized codes rank 0 so they surface.
flag_priority <- c("E", "M", "ZM", "B", "Z", "L", "D", "Dry", "ND", "none")

q_daily <- q_sub %>%
  filter(!is.na(datetime), year(datetime) >= 2021, year(datetime) <= 2024) %>%
  distinct(region, siteId, datetime, .keep_all = TRUE) %>%
  mutate(
    date      = as_date(datetime),
    flag      = coalesce(QAQC, na_flag, "none"),
    flag_rank = coalesce(match(flag, flag_priority), 0L)
  ) %>%
  group_by(region, siteId, date) %>%
  summarise(
    n_rows    = n(),
    n_obs     = sum(!is.na(Q_Ls)),
    Q_max     = if (n_obs > 0) max(Q_Ls, na.rm = TRUE) else NA_real_,   # before Q_Ls is overwritten by the mean
    Q_Ls      = if (n_obs > 0) mean(Q_Ls, na.rm = TRUE) else NA_real_,
    n_dry     = sum(wetdry %in% "dry"),
    wetdry    = case_when(
      all(is.na(wetdry)) ~ NA_character_,
      n_dry == 0         ~ "wet",
      n_dry == n_rows    ~ "dry",
      TRUE               ~ "mixed"
    ),
    n_flagged = sum(flag != "none"),
    flag_plot = flag[which.min(flag_rank)],
    flag_all  = paste(sort(unique(flag)), collapse = ";"),
    .groups = "drop"
  )

# Check: days per site and date range
q_daily %>%
  group_by(region, siteId) %>%
  summarise(n_days = n(), start = min(date), end = max(date),
            n_days_with_q = sum(!is.na(Q_Ls)),
            n_dry_days = sum(wetdry %in% "dry"), .groups = "drop")

# Check: flag codes that are not in the priority list
q_daily %>%
  filter(!flag_plot %in% flag_priority) %>%
  count(region, flag_plot)

# Check: max daily mean discharge by site
q_daily %>%
  group_by(region) %>%
  summarise(max_daily_mean = max(Q_Ls, na.rm = TRUE),
            max_daily_max  = max(Q_max, na.rm = TRUE), .groups = "drop")

# 4.0 Hydrographs colored by flag ----------------------------------------------
flag_levels <- c("none", "Z", "M", "ZM", "B", "L", "D", "E", "Dry", "ND")

flag_cols <- c(
  "none" = "steelblue",
  "Z"    = "#E69F00",   # outlier, z > 3
  "M"    = "#D55E00",   # modeled > 2x max measured
  "ZM"   = "#CC79A7",   # outlier and modeled > 2x max measured
  "B"    = "#F0E442",   # above bankfull
  "L"    = "#009E73",   # modeled < half min measured
  "D"    = "grey40",    # within 2 cm sensor error (dry)
  "E"    = "black",     # sensor error code
  "Dry"  = "grey40",    # GBJ, JDR dry
  "ND"   = "#999999"    # GBJ, JDR no data
)

site_labels <- c(
  GBJ = "GBJ (GSS01)", JDR = "JDR (JSS01)", KNZ = "KNZ (SFM01)", PRF = "PRF (PRM01)",
  TAL = "TAL (TLM01)", WHR = "WHR (WHM01)", YMR = "YMR (ENM01)"
)

q_plot <- q_daily %>%
  filter(!is.na(Q_Ls), Q_Ls > 0) %>%      # log scale needs Q > 0
  mutate(
    flag_plot = factor(flag_plot, levels = flag_levels),
    region    = factor(region, levels = names(site_labels), labels = site_labels)
  ) %>%
  arrange(flag_plot != "none", flag_plot)  # draw flagged days on top

p_hydro <- ggplot(q_plot, aes(x = date, y = Q_Ls)) +
  geom_line(color = "grey80", lwd = 0.3) +
  geom_point(aes(color = flag_plot), size = 1.2, alpha = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = flag_cols, drop = TRUE, name = "Flag") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  facet_wrap(~region, ncol = 2, scales = "free") +
  labs(x = NULL, y = "Daily mean discharge (L/s)") +
  theme_bw() +
  theme(
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 10),
    strip.text      = element_text(size = 11, face = "bold"),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    legend.position = "right"
  )

p_hydro

dir.create("figures", showWarnings = FALSE)
ggsave("figures/flow_hydrographs_by_flag.png", p_hydro,
       width = 11, height = 9, dpi = 300)

# Count of days by flag, and days hidden by the log scale (zero or missing flow)
q_daily %>%
  group_by(region, flag_plot) %>%
  summarise(n_days = n(),
            n_plotted = sum(!is.na(Q_Ls) & Q_Ls > 0),
            n_hidden = n_days - n_plotted, .groups = "drop") %>%
  print(n = Inf)

# 7.0 Chemistry sample days on the flow duration curves ------------------------
chem <- read_csv("data/allSites_chem_Q_data.csv", show_col_types = FALSE) %>%
  filter(is_outlet) %>%
  transmute(
    siteId, watershed, date = as_date(date),
    Q_chem = as.numeric(Q_Ls),
    Q_flag, Q_partial_flag
  ) %>%
  distinct(siteId, date, .keep_all = TRUE)    # one row per sampling day

# 7.1 Flow duration curves, ties ranked by their average
fdc <- q_daily %>%
  filter(!is.na(Q_Ls)) %>%
  group_by(region, siteId) %>%
  mutate(exceed_pct = 100 * rank(-Q_Ls, ties.method = "average") / (n() + 1)) %>%
  ungroup()

# 7.2 Attach exceedance to each sample day
samples <- chem %>%
  inner_join(fdc %>% select(siteId, date, region, Q_daily = Q_Ls, exceed_pct),
             by = c("siteId", "date")) %>%
  mutate(flagged = if_else(is.na(Q_flag), "no flag", "flagged"))

# Check: samples lost in the join, by site
chem %>%
  anti_join(fdc, by = c("siteId", "date")) %>%
  count(siteId, watershed, name = "n_no_daily_flow")

# Check: chem Q_Ls against our daily mean (differences point to a date offset or M removal)
samples %>%
  mutate(ratio = Q_chem / Q_daily) %>%
  group_by(siteId) %>%
  summarise(n = n(), median_ratio = median(ratio, na.rm = TRUE),
            n_off_by_10pct = sum(abs(ratio - 1) > 0.1, na.rm = TRUE),
            .groups = "drop")

# 7.3 Where did samples fall? Exceedance quantiles by site
samples %>%
  group_by(region, siteId) %>%
  summarise(
    n_samples   = n(),
    min_exceed  = min(exceed_pct),
    q25         = quantile(exceed_pct, 0.25),
    median      = median(exceed_pct),
    q75         = quantile(exceed_pct, 0.75),
    max_exceed  = max(exceed_pct),
    n_zero_flow = sum(Q_daily == 0),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# 7.4 Figure: curve with a vertical line at each sample's exceedance
site_labels <- c(
  GBJ = "GBJ (GSS01)", JDR = "JDR (JSS01)", KNZ = "KNZ (SFM01)", PRF = "PRF (PRM01)",
  TAL = "TAL (TLM01)", WHR = "WHR (WHM01)", YMR = "YMR (ENM01)"
)

fdc_plot <- fdc %>%
  filter(Q_Ls > 0) %>%
  mutate(region = factor(region, levels = names(site_labels), labels = site_labels))

samp_plot <- samples %>%
  mutate(region = factor(region, levels = names(site_labels), labels = site_labels))

p_fdc <- ggplot() +
  geom_vline(data = samp_plot,
             aes(xintercept = exceed_pct, color = flagged),
             alpha = 0.5, lwd = 0.4) +
  geom_line(data = fdc_plot, aes(x = exceed_pct, y = Q_Ls),
            color = "steelblue", lwd = 0.75) +
  scale_y_log10() +
  scale_color_manual(values = c("no flag" = "grey30", "flagged" = "#D55E00"),
                     name = "Flow flag on\nsample day") +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 1.2))) +
  facet_wrap(~region, ncol = 2, scales = "free_y") +
  labs(x = "Percent of days flow is exceeded", y = "Daily mean discharge (L/s)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"))

p_fdc

ggsave("figures/fdc_with_chem_samples.png", p_fdc, width = 11, height = 9, dpi = 300)

# 7.5 Accounting: where did days go, by site?
q_daily %>%
  group_by(region, siteId) %>%
  summarise(
    n_days_total   = n(),
    n_days_no_data = sum(is.na(Q_Ls)),
    n_days_zero    = sum(Q_Ls == 0, na.rm = TRUE),
    n_days_pos     = sum(Q_Ls > 0, na.rm = TRUE),
    pct_zero_of_data = round(100 * n_days_zero / (n_days_zero + n_days_pos), 1),
    .groups = "drop"
  ) %>%
  print(width = Inf)

# 7.6 Same figure, with the zero-flow zone shaded
zero_zone <- fdc %>%
  group_by(region) %>%
  summarise(
    xmin = 100 * sum(Q_Ls > 0) / (n() + 1),
    pct_zero = round(100 * mean(Q_Ls == 0), 0),
    .groups = "drop"
  ) %>%
  filter(pct_zero > 0) %>%
  mutate(region = factor(region, levels = names(site_labels), labels = site_labels),
         xmax = 100)

p_fdc2 <- ggplot() +
  geom_rect(data = zero_zone,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.6) +
  geom_text(data = zero_zone,
            aes(x = (xmin + xmax) / 2, y = Inf, label = paste0("no flow\n", pct_zero, "%")),
            vjust = 1.3, size = 3.2, color = "grey30") +
  geom_vline(data = samp_plot,
             aes(xintercept = exceed_pct, color = flagged),
             alpha = 0.5, lwd = 0.4) +
  geom_line(data = fdc_plot, aes(x = exceed_pct, y = Q_Ls),
            color = "steelblue", lwd = 0.75) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("no flag" = "grey30", "flagged" = "#D55E00"),
                     name = "Flow flag on\nsample day") +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 1.2))) +
  facet_wrap(~region, ncol = 2, scales = "free_y") +
  labs(x = "Percent of days flow is exceeded", y = "Daily mean discharge (L/s)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"))

p_fdc2

ggsave("figures/fdc_with_chem_samples_zero_shaded.png", p_fdc2,
       width = 11, height = 9, dpi = 300)

# 8.0 Flow duration curves, all sites on one plot --------------------------------
site_cols <- c(
  "GBJ (GSS01)" = "#0072B2", "JDR (JSS01)" = "#56B4E9", "KNZ (SFM01)" = "#D55E00",
  "PRF (PRM01)" = "#009E73", "TAL (TLM01)" = "#CC79A7", "WHR (WHM01)" = "#E69F00",
  "YMR (ENM01)" = "#000000"
)

fdc_all <- fdc %>%
  filter(Q_Ls > 0) %>%
  mutate(site = factor(region, levels = names(site_labels), labels = site_labels))

# Raw discharge
p_all <- ggplot(fdc_all, aes(x = exceed_pct, y = Q_Ls, color = site)) +
  geom_line(lwd = 0.75) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = site_cols, name = NULL) +
  labs(x = "Percent of days flow is exceeded", y = "Daily mean discharge (L/s)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        legend.text = element_text(size = 11))

p_all

# Same curves, discharge divided by each site's median flow on days with flow
fdc_norm <- fdc_all %>%
  group_by(site) %>%
  mutate(Q_rel = Q_Ls / median(Q_Ls)) %>%
  ungroup()

p_norm <- ggplot(fdc_norm, aes(x = exceed_pct, y = Q_rel, color = site)) +
  geom_line(lwd = 0.75) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = site_cols, name = NULL) +
  labs(x = "Percent of days flow is exceeded",
       y = "Discharge / median discharge on flowing days") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        legend.text = element_text(size = 11))

p_norm

ggsave("figures/fdc_all_sites_raw.png",  p_all,  width = 8, height = 6, dpi = 300)
ggsave("figures/fdc_all_sites_norm.png", p_norm, width = 8, height = 6, dpi = 300)

# 9.0 Flow duration curves from all observations, with chem sample days --------
fdc_obs <- q_sub %>%
  filter(!is.na(datetime), !is.na(Q_Ls)) %>%
  distinct(region, siteId, datetime, .keep_all = TRUE) %>%
  select(region, siteId, Q_Ls)

# Rank every observation, then thin to 1,500 points per site for plotting only
fdc_obs_plot <- fdc_obs %>%
  group_by(region, siteId) %>%
  arrange(desc(Q_Ls), .by_group = TRUE) %>%
  mutate(exceed_pct = 100 * rank(-Q_Ls, ties.method = "average") / (n() + 1)) %>%
  slice(unique(round(seq(1, n(), length.out = 1500)))) %>%
  ungroup() %>%
  filter(Q_Ls > 0) %>%
  mutate(region = factor(region, levels = names(site_labels), labels = site_labels))

# Exceedance of each sample day's discharge within the all-observation distribution
# Ties (zero flow) get the midpoint of the tied block
obs_by_site <- split(fdc_obs$Q_Ls, fdc_obs$siteId)

samples_obs <- samples %>%
  mutate(
    exceed_obs = map2_dbl(siteId, Q_daily, function(s, q) {
      x <- obs_by_site[[s]]
      100 * (mean(x > q) + 0.5 * mean(x == q))
    }),
    region = factor(region, levels = names(site_labels), labels = site_labels)
  )

# Zero-flow zone from all observations
zero_zone_obs <- fdc_obs %>%
  group_by(region) %>%
  summarise(xmin = 100 * mean(Q_Ls > 0), pct_zero = round(100 * mean(Q_Ls == 0), 0),
            .groups = "drop") %>%
  filter(pct_zero > 0) %>%
  mutate(region = factor(region, levels = names(site_labels), labels = site_labels),
         xmax = 100)

p_fdc_obs <- ggplot() +
  geom_rect(data = zero_zone_obs,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.6) +
  geom_text(data = zero_zone_obs,
            aes(x = (xmin + xmax) / 2, y = Inf, label = paste0("no flow\n", pct_zero, "%")),
            vjust = 1.3, size = 3.2, color = "grey30") +
  geom_vline(data = samples_obs,
             aes(xintercept = exceed_obs, color = flagged),
             alpha = 0.5, lwd = 0.4) +
  geom_line(data = fdc_obs_plot, aes(x = exceed_pct, y = Q_Ls),
            color = "steelblue", lwd = 0.75) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("no flag" = "grey30", "flagged" = "#D55E00"),
                     name = "Flow flag on\nsample day") +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 1.2))) +
  facet_wrap(~region, ncol = 2, scales = "free_y") +
  labs(x = "Percent of time flow is exceeded", y = "Discharge (L/s)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold"))

p_fdc_obs

ggsave("figures/fdc_all_obs_with_chem_samples.png", p_fdc_obs,
       width = 11, height = 9, dpi = 300)

# Compare sample exceedance: daily-curve method versus all-observation method
samples_obs %>%
  group_by(region) %>%
  summarise(median_exceed_daily = median(exceed_pct),
            median_exceed_obs   = median(exceed_obs),
            max_abs_diff        = max(abs(exceed_pct - exceed_obs)),
            .groups = "drop")

# 10.0 Chemistry versus flow, log-log regression, all outlet sites -------------
# Needs fdc (7.1), site_labels (7.4), and site_cols (8.0) in the session

# Switches
x_var           <- "Q_rel"   # "Q_rel" = discharge / site median on flowing days
# "flow_quantile" = percentile (0 = lowest flow, 100 = highest)
drop_bad_Q_days <- TRUE      # TRUE drops sample days whose chem-file Q_flag has M, B, or E
# from the regressions (they stay on the plots)
min_n_per_site  <- 8         # minimum flowing samples to fit a site regression

x_lab <- if (x_var == "Q_rel") {
  "Discharge / site median discharge on flowing days"
} else {
  "Flow quantile (0 = lowest flow, 100 = highest)"
}
slope_lab <- if (x_var == "Q_rel") {
  "log-log slope (0 = chemostatic, < 0 dilution, > 0 enrichment)"
} else {
  "log-log slope against flow quantile"
}

solutes <- c("F_avg", "Cl_avg", "NO3_avg", "Br_avg", "SO4_avg", "Na_avg", "Ca_avg",
             "B_avg", "Mg_avg", "Si_avg", "K_avg", "Sr_avg",
             "SRPugL", "NH4ugL", "NO3ugL", "NPOCmgL")
optical <- c("BIX", "HIX", "FI", "A254", "A300", "e2e3", "s275_295", "s350_400")

# 10.1 Join chemistry to each sample day's flow
site_Qmed <- fdc %>%
  filter(Q_Ls > 0) %>%
  group_by(siteId) %>%
  summarise(Q_med = median(Q_Ls), .groups = "drop")

chem_q <- read_csv("data/allSites_chem_Q_data.csv", show_col_types = FALSE) %>%
  filter(is_outlet) %>%
  rename(aims_region = region) %>%    # chem file has its own region column (GP, MW, SE)
  mutate(date = as_date(date)) %>%
  distinct(siteId, date, .keep_all = TRUE) %>%
  inner_join(fdc %>% select(siteId, date, region, Q_daily = Q_Ls, exceed_pct),
             by = c("siteId", "date")) %>%
  left_join(site_Qmed, by = "siteId") %>%
  mutate(
    flow_quantile = 100 - exceed_pct,
    Q_rel         = Q_daily / Q_med,
    zero_flow     = Q_daily == 0,
    bad_Q         = str_detect(replace_na(as.character(Q_flag), ""), "M|B|E"),
    region = factor(region, levels = names(site_labels), labels = site_labels)
  )

# Check: samples per site, zero-flow samples, and sample days with suspect Q
chem_q %>%
  group_by(region) %>%
  summarise(n_samples = n(), n_zero_flow = sum(zero_flow), n_bad_Q = sum(bad_Q),
            .groups = "drop") %>%
  print(n = Inf)

# 10.2 Long format, one row per sample and analyte
chem_long <- chem_q %>%
  select(siteId, region, date, flow_quantile, Q_rel, zero_flow, bad_Q,
         all_of(c(solutes, optical))) %>%
  pivot_longer(all_of(c(solutes, optical)), names_to = "analyte", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value)) %>%
  group_by(siteId, analyte) %>%
  mutate(value_rel = value / median(value[value > 0])) %>%   # concentration relative to site median
  ungroup()

# Check: samples per site, and values that log axes will hide
chem_long %>% count(region) %>% print(n = Inf)
chem_long %>% group_by(analyte) %>%
  summarise(n = n(), n_zero_or_neg = sum(value <= 0), .groups = "drop") %>%
  print(n = Inf)

# Check: extreme values relative to the site median (look for single bad samples)
chem_long %>%
  filter(value_rel > 8 | value_rel < 0.1) %>%
  arrange(siteId, date, analyte) %>%
  select(siteId, date, flow_quantile, analyte, value, value_rel) %>%
  print(n = Inf)

# 10.3 Data for plots and regressions
# Points: every sample with a positive x. Fits: flowing days, and good Q if the switch is on.
chem_cq <- chem_long %>%
  mutate(x = .data[[x_var]]) %>%
  filter(!is.na(x), x > 0, value > 0)

chem_fit <- chem_cq %>%
  filter(!zero_flow, !(drop_bad_Q_days & bad_Q))

# 10.4 Plot function: points, one thin regression line per site, one pooled line
plot_chem_cq <- function(analytes, yvar = "value_rel", ylab = "Concentration / site median") {
  keep <- function(d) {
    d %>%
      filter(analyte %in% analytes, !is.na(.data[[yvar]]), .data[[yvar]] > 0) %>%
      mutate(analyte = factor(analyte, levels = analytes))
  }
  pts <- keep(chem_cq)
  fit <- keep(chem_fit)
  site_fit <- fit %>%
    group_by(analyte, region) %>%
    filter(n() >= min_n_per_site) %>%
    ungroup()
  
  ggplot(pts, aes(x = x, y = .data[[yvar]])) +
    geom_point(aes(color = region, shape = zero_flow), size = 1.4, alpha = 0.7) +
    geom_smooth(data = site_fit, aes(color = region, group = region),
                method = "lm", formula = y ~ x, se = FALSE, lwd = 0.4) +
    geom_smooth(data = fit, method = "lm", formula = y ~ x, se = FALSE,
                color = "black", lwd = 1) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = site_cols, name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1),
                       labels = c("FALSE" = "flowing", "TRUE" = "no flow"), name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, shape = 16))) +
    facet_wrap(~analyte, scales = "free_y", ncol = 4) +
    labs(x = x_lab, y = ylab) +
    theme_bw() +
    theme(axis.title      = element_text(size = 14),
          axis.text       = element_text(size = 10),
          strip.text      = element_text(size = 11, face = "bold"),
          legend.text     = element_text(size = 11),
          legend.position = "bottom")
}

# 10.5 Figures
p_solutes_rel <- plot_chem_cq(solutes, "value_rel", "Concentration / site median")
p_solutes_raw <- plot_chem_cq(solutes, "value",     "Concentration (units as reported)")
p_optical_rel <- plot_chem_cq(optical, "value_rel", "Value / site median")

p_solutes_rel
p_solutes_raw
p_optical_rel

dir.create("figures", showWarnings = FALSE)
ggsave("figures/chem_loglog_solutes_rel.png", p_solutes_rel, width = 13, height = 11, dpi = 300)
ggsave("figures/chem_loglog_solutes_raw.png", p_solutes_raw, width = 13, height = 11, dpi = 300)
ggsave("figures/chem_loglog_optical_rel.png", p_optical_rel, width = 13, height = 7,  dpi = 300)

# 10.6 Regression slopes in log-log space, by analyte and site, plus pooled
cq_fit <- function(d, y) {
  d <- d %>% mutate(ly = log10(.data[[y]]), lx = log10(x))
  m <- lm(ly ~ lx, data = d)
  broom::tidy(m, conf.int = TRUE) %>%
    filter(term == "lx") %>%
    transmute(n = nobs(m), slope = estimate, ci_low = conf.low, ci_high = conf.high,
              p = p.value, r2 = summary(m)$r.squared)
}

cq_site <- chem_fit %>%
  group_by(analyte, region) %>%
  filter(n() >= min_n_per_site) %>%
  group_modify(~cq_fit(.x, "value")) %>%
  ungroup() %>%
  mutate(p_BH = p.adjust(p, method = "BH"))     # adjusts across all site-by-analyte tests

cq_pooled <- chem_fit %>%
  filter(!is.na(value_rel), value_rel > 0) %>%
  group_by(analyte) %>%
  group_modify(~cq_fit(.x, "value_rel")) %>%
  ungroup() %>%
  mutate(region = factor("Pooled"), p_BH = NA_real_)

cq_all <- bind_rows(cq_site %>% mutate(region = as.character(region)),
                    cq_pooled %>% mutate(region = as.character(region))) %>%
  mutate(
    region  = factor(region, levels = c(names(site_cols), "Pooled")),
    analyte = factor(analyte, levels = rev(c(solutes, optical)))
  )

cq_all %>%
  arrange(analyte, region) %>%
  mutate(across(c(slope, ci_low, ci_high, r2), ~round(.x, 2)),
         across(c(p, p_BH), ~signif(.x, 2))) %>%
  print(n = Inf, width = Inf)

# 10.7 Slope summary figure
p_slopes <- ggplot(cq_all, aes(x = slope, y = analyte, color = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high),
                  position = position_dodge(width = 0.7), size = 0.3) +
  scale_color_manual(values = c(site_cols, "Pooled" = "black"), name = NULL) +
  labs(x = slope_lab, y = NULL) +
  theme_bw() +
  theme(axis.title      = element_text(size = 14),
        axis.text       = element_text(size = 10),
        legend.text     = element_text(size = 11),
        legend.position = "bottom")

p_slopes

ggsave("figures/chem_loglog_slopes.png", p_slopes, width = 9, height = 11, dpi = 300)