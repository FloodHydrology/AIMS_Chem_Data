#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Stream Elasticity Quantification
# Date: 10/27/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Explore elasticity across sites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Approach
#  Step 1: Download all flow data
#  Step 2: Estimate daily flow and ASDN metrics
#  Step 3: Explore B scaling factor
#  Step 4: Explore variation in flow vs variation in ASDN
#  Step 5: Make badass data viz

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
    Q_Ls      = if (all(is.na(Q_Ls))) NA_real_ else mean(Q_Ls, na.rm = TRUE),
    flag_plot = flag[which.min(flag_rank)],
    flag_all  = paste(sort(unique(flag)), collapse = ";"),
    .groups = "drop"
  )

# Check: flag codes that are not in the priority list
q_daily %>%
  filter(!flag_plot %in% flag_priority) %>%
  count(region, flag_plot)

# 4.0 Load ASDN and join to daily Q --------------------------------------------
# TODO: replace with your actual ASDN read. Needs columns:
# date, watershed, percent_network_wet, active_network_length_m
asdn <- read_csv("data/streamDAG_metrics.csv") |> select(watershed, date, percent_network_wet, active_network_length_m)

# ASDN as a fraction of the full network (length at 100% wet)
asdn <- asdn %>%
  group_by(watershed) %>%
  mutate(network_len_m = max(active_network_length_m, na.rm = TRUE),
         asdn_frac     = active_network_length_m / network_len_m) %>%
  ungroup()

q_asdn <- q_daily %>%
  left_join(asdn, by = c("region" = "watershed", "date"))

# Check: paired days by site
q_asdn %>%
  group_by(region, siteId) %>%
  summarise(n_days    = n(),
            n_no_asdn = sum(is.na(active_network_length_m)),
            n_no_q    = sum(is.na(Q_Ls)),
            n_both    = sum(!is.na(Q_Ls) & !is.na(active_network_length_m)),
            n_zero_q  = sum(Q_Ls == 0, na.rm = TRUE),
            .groups = "drop")

# 5.0 Convert discharge to specific discharge (mm/day) -------------------------
wshed <- read_csv("data/aims_wshed_characteristics.csv")

# Area codes differ from the Q and ASDN codes for four sites. Match by site name.
area_key <- wshed %>%
  transmute(
    region = case_when(
      watershedme == "Shane" ~ "SHN",
      TRUE ~ recode(watershed, KMZ = "KNZ", PRA = "PRF", DCR = "DRC")
    ),
    area_km2
  )

# Checks: all Q sites have an area, and no duplicate keys
anti_join(distinct(q_daily, region), area_key, by = "region")
stopifnot(!anyDuplicated(area_key$region))

# 1 L/s = 86.4 m3/day. 1 km2 = 1e6 m2. So mm/day = L/s * 0.0864 / km2.
q_asdn <- q_asdn %>%
  left_join(area_key, by = "region") %>%
  mutate(Q_mm_d = Q_Ls * 0.0864 / area_km2)

# Check: specific discharge by site
q_asdn %>%
  group_by(region) %>%
  summarise(area_km2   = first(area_km2),
            Q_mm_d_p05 = quantile(Q_mm_d, 0.05, na.rm = TRUE),
            Q_mm_d_med = median(Q_mm_d, na.rm = TRUE),
            Q_mm_d_p95 = quantile(Q_mm_d, 0.95, na.rm = TRUE),
            .groups = "drop")

# 6.0 Plot Q vs ASDN in log-log space ------------------------------------------
# Zero-flow days plot at half the smallest nonzero Q per site. Display only.
# Do not fit on these values.
q_plot <- q_asdn %>%
  filter(!is.na(Q_mm_d), !is.na(asdn_frac), asdn_frac > 0) %>%
  group_by(region) %>%
  mutate(
    q_floor   = min(Q_mm_d[Q_mm_d > 0]) / 2,
    zero_q    = Q_mm_d == 0,
    Q_plot    = if_else(zero_q, q_floor, Q_mm_d),
    saturated = percent_network_wet == 100
  ) %>%
  ungroup()

q_plot %>%
  ggplot(aes(x = Q_plot, y = asdn_frac)) +
  geom_point(aes(color = saturated, shape = zero_q), alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkorange"),
                     labels = c("FALSE" = "Some STICs dry", "TRUE" = "All STICs wet"),
                     name   = NULL) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 2),
                     labels = c("FALSE" = "Flowing", "TRUE" = "Zero flow (plotted at floor)"),
                     name   = NULL) +
  facet_wrap(~region, scales = "free_x") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Daily mean discharge (mm/day)",
       y = "ASDN (fraction of network)") +
  theme_bw() +
  theme(axis.title      = element_text(size = 14),
        axis.text       = element_text(size = 10),
        legend.position = "bottom")
