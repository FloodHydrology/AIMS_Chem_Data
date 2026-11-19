#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Budyko analysis of the AIMS watersheds
# Date: 11/18/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Place each watershed-water year in Budyko space (evaporative index vs.
#          aridity index), with point size showing how wet the stream network
#          was over that water year (mean daily ASDN).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Approach
#  - Water year (WY): 1 Oct - 30 Sep, named by the ending year.
#  - P and PET: gridMET daily precipitation and grass reference ET (Penman-
#    Monteith; Abatzoglou 2013, ~4 km) at each watershed outlet, from the public
#    THREDDS server (point queries). Cached to output/ after the first download.
#  - Q: annual runoff (mm) = mean daily Q_mm_day x days in the WY, from the merged
#    streamDAG/discharge file (E. Seybold). Assumes missing days resemble observed
#    days (see coverage notes below).
#  - Aridity index = PET / P. Evaporative index = AET / P = 1 - Q / P (annual
#    storage change assumed ~0).
#  - Network wetness = mean daily ASDN_pct over the WY (% of maximum active length).
#  - Reference curve: Budyko (1974), E/P = sqrt(phi * tanh(1/phi) * (1 - exp(-phi))),
#    with energy (E = PET) and water (E = P) limits.
#
# Inputs
#  data/AIMS_Q_connectivity_merge/output/streamDAG_Q_data.csv  daily Q and ASDN
#  data/AIMS_Q_connectivity_merge/aims_data/*_Q/*ReadME*.xlsx  outlet coordinates
#  data/site_colors.xlsx                                       site colors
#  aims_data/<SE site>_METS/METS_*.csv                         on-site rain gauges
#     (PRF, TAL, WHR only; used to check gridMET P)
#
# Notes and caveats
#  - Q coverage: WYs with >= 85% of days with Q are "complete" (filled points);
#    50-85% are "incomplete" (open points); < 50% are dropped. Missing Q days are
#    NOT dry days in general: GBJ WY2023 lacks most of the May-June snowmelt peak
#    (likely readings above the rating curve removed in the merge), and YMR gaps
#    are mostly winter. Incomplete years, and GBJ WY2023, likely underestimate Q
#    (overestimate the evaporative index).
#  - GBJ area is unconfirmed (16.2 vs. 14.83 km2); with 14.83 km2 GBJ runoff is
#    ~9% higher.
#  - KNZ and YMR outlet coordinates are not in the local ReadMes; taken from the
#    AIMS HydroShare resources (KNZ SFM01 39.092281, -96.58719; YMR outlet
#    37.56442, -96.49106).
#  - gridMET is a ~4 km grid; small watersheds sit within one or two cells.
#    Check against on-site gauges (section 2.1): WY totals within 1% at TAL
#    (WY2024) and WHR (WY2022).
#  - PRF WY2022 runoff is 84% of P (evaporative index 0.16), implausible for a
#    forested Southeast watershed: likely a discharge (rating curve) or area
#    issue. Question for E. Seybold, along with the PRF 2021 ASDN band.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1.0 Setup workspace ----------------------------------------------------------
#Clear workspace
remove(list = ls())

#Load libraries of interest
library(tidyverse)
library(lubridate)
library(readxl)
library(ggrepel)

#Define file paths
merge_file   <- "data/AIMS_Q_connectivity_merge/output/streamDAG_Q_data.csv"
q_dir        <- "data/AIMS_Q_connectivity_merge/aims_data"
colors_file  <- "data/site_colors.xlsx"
fig_dir      <- "figures"
tab_dir      <- "output"
gridmet_file <- file.path(tab_dir, "climate_daily_outlets.csv")   # P and PET cache

#Analysis settings
wy_start       <- as_date("2020-10-01")   # first WY with any data: 2021
wy_end         <- as_date("2025-09-30")
cov_complete   <- 0.85   # Q coverage for a "complete" WY
cov_min        <- 0.50   # minimum Q coverage to show a WY at all
asdn_cov_min   <- 0.50   # minimum ASDN coverage for the wetness value
refresh_gridmet <- FALSE # TRUE to re-download gridMET

water_year <- function(d) year(d) + (month(d) >= 10)

# 1.1 Sites, order, colors -----------------------------------------------------
site_info <- tribble(
  ~watershed, ~region,
  "KNZ", "Great Plains", "YMR", "Great Plains",
  "GBJ", "Mountain West", "JDR", "Mountain West",
  "PRF", "Southeast", "TAL", "Southeast", "WHR", "Southeast"
)
site_order <- site_info$watershed
col_line <- read_excel(colors_file) %>% select(site, hex2) %>% deframe()
col_line <- col_line[site_order]

# 1.2 Outlet coordinates -------------------------------------------------------
# From the Site_Info sheet of each discharge ReadMe where present
readme_coords <- map_dfr(site_order, function(ws) {
  f <- list.files(file.path(q_dir, paste0(ws, "_Q")), pattern = "ReadME.*\\.xlsx$",
                  full.names = TRUE)
  if (length(f) == 0 || !"Site_Info" %in% excel_sheets(f[1])) return(NULL)
  read_excel(f[1], sheet = "Site_Info") %>%
    transmute(watershed = ws, siteId, latitude, longitude, coord_source = "discharge ReadMe")
})

# KNZ and YMR: from the AIMS HydroShare resources (see header)
outlets <- bind_rows(
  readme_coords,
  tribble(~watershed, ~siteId, ~latitude, ~longitude, ~coord_source,
          "KNZ", "SFM01", 39.092281, -96.58719, "HydroShare (AIMS KNZ discharge)",
          "YMR", "ENM01", 37.56442,  -96.49106, "HydroShare (AIMS YMR outlet)")
) %>%
  filter(watershed %in% site_order) %>%
  distinct(watershed, .keep_all = TRUE) %>%
  arrange(factor(watershed, levels = site_order))
stopifnot(setequal(outlets$watershed, site_order))
print(outlets)

# 2.0 gridMET P and PET --------------------------------------------------------
# Primary source: gridMET per-year files (pr_YYYY.nc, pet_YYYY.nc) through the
# THREDDS NetCDF Subset Service. The CSV output is NOT unpacked, so each value is
# multiplied by the scale_factor/add_offset read from the file's own metadata.
# The server is intermittently unavailable (HTTP 503), so every request is retried
# and each site-variable-year is cached in output/gridmet_cache/ (reruns resume).
# Fallback: if gridMET cannot be completed for every site, ALL sites switch to
# Open-Meteo's historical archive (ERA5-based precipitation and FAO-56 grass
# reference ET, ~9-25 km) so that sources are never mixed.
thredds   <- "https://thredds.northwestknowledge.net/thredds"
cache_dir <- file.path(tab_dir, "gridmet_cache")
dir.create(cache_dir, showWarnings = FALSE)

fetch_text <- function(url, ok, tries = 20, wait = 30) {
  for (i in seq_len(tries)) {
    txt <- tryCatch(suppressWarnings(readLines(url(url), warn = FALSE)),
                    error = function(e) character(0))
    if (length(txt) && ok(txt)) return(txt)
    if (i < tries) Sys.sleep(wait)
  }
  NULL
}

# Data variable name and packing for one gridMET per-year file
gridmet_meta <- function(v, yr) {
  das <- fetch_text(sprintf("%s/dodsC/MET/%s/%s_%d.nc.das", thredds, v, v, yr),
                    ok = function(t) any(grepl("scale_factor", t)))
  if (is.null(das)) return(NULL)
  starts <- grep("^    [A-Za-z_]+ \\{", das)
  names_ <- sub("^    ([A-Za-z_]+) \\{.*", "\\1", das[starts])
  k      <- which(!names_ %in% c("lat", "lon", "day", "crs", "NC_GLOBAL"))[1]
  block  <- das[starts[k]:(c(starts, length(das))[k + 1])]
  num    <- function(key) {
    x <- grep(key, block, value = TRUE)
    if (length(x)) as.numeric(sub(".* ([-0-9.eE]+);", "\\1", x[1])) else NA_real_
  }
  list(var = names_[k], scale = coalesce(num("scale_factor"), 1),
       offset = coalesce(num("add_offset"), 0))
}

gridmet_year <- function(ws, lat, lon, v, yr) {
  cache <- file.path(cache_dir, sprintf("%s_%s_%d.csv", ws, v, yr))
  if (file.exists(cache)) return(read_csv(cache, show_col_types = FALSE))
  meta <- gridmet_meta(v, yr)
  if (is.null(meta)) return(NULL)
  t0 <- max(wy_start, as_date(sprintf("%d-01-01", yr)))
  t1 <- min(wy_end,   as_date(sprintf("%d-12-31", yr)))
  url <- sprintf("%s/ncss/MET/%s/%s_%d.nc?var=%s&latitude=%s&longitude=%s&time_start=%sT00:00:00Z&time_end=%sT00:00:00Z&accept=csv",
                 thredds, v, v, yr, meta$var, lat, lon, t0, t1)
  n_days <- as.numeric(t1 - t0) + 1
  txt <- fetch_text(url, ok = function(t) grepl("^time,", t[1]) && length(t) >= n_days + 1)
  if (is.null(txt)) return(NULL)
  d <- read_csv(I(txt), show_col_types = FALSE) %>%
    transmute(date = as_date(time), value = .[[4]] * meta$scale + meta$offset)
  write_csv(d, cache)
  d
}

get_gridmet <- function() {
  yrs <- year(wy_start):year(wy_end)
  out <- pmap(outlets, function(watershed, latitude, longitude, ...) {
    message("gridMET: ", watershed)
    vals <- map(c(P_mm = "pr", PET_mm = "pet"), function(v) {
      parts <- map(yrs, ~gridmet_year(watershed, latitude, longitude, v, .x))
      if (any(map_lgl(parts, is.null))) NULL else bind_rows(parts)
    })
    if (any(map_lgl(vals, is.null))) return(NULL)
    full_join(rename(vals$P_mm, P_mm = value), rename(vals$PET_mm, PET_mm = value),
              by = "date") %>%
      mutate(watershed = watershed, .before = 1)
  })
  if (any(map_lgl(out, is.null))) NULL else bind_rows(out) %>% mutate(source = "gridMET")
}

get_openmeteo <- function() {
  pmap_dfr(outlets, function(watershed, latitude, longitude, ...) {
    message("Open-Meteo: ", watershed)
    url <- sprintf("https://archive-api.open-meteo.com/v1/archive?latitude=%s&longitude=%s&start_date=%s&end_date=%s&daily=precipitation_sum,et0_fao_evapotranspiration&timezone=auto&format=csv",
                   latitude, longitude, wy_start, wy_end)
    txt <- fetch_text(url, ok = function(t) any(grepl("^time,", t)), tries = 5, wait = 20)
    if (is.null(txt)) stop("Open-Meteo download failed for ", watershed)
    read_csv(I(txt[grep("^time,", txt):length(txt)]), show_col_types = FALSE) %>%
      transmute(watershed, date = as_date(time), P_mm = .[[2]], PET_mm = .[[3]],
                source = "Open-Meteo (ERA5)")
  })
}

if (!file.exists(gridmet_file) || refresh_gridmet) {
  climate_daily <- get_gridmet()
  if (is.null(climate_daily)) {
    message("gridMET incomplete; using Open-Meteo for all sites")
    climate_daily <- get_openmeteo()
  }
  write_csv(climate_daily, gridmet_file)
} else {
  climate_daily <- read_csv(gridmet_file, show_col_types = FALSE)
}
climate_source <- unique(climate_daily$source)
cat("\nClimate source:", climate_source, "\n")

climate_wy <- climate_daily %>%
  mutate(WY = water_year(date)) %>%
  group_by(watershed, WY) %>%
  summarise(n_clim = n(), P_mm = sum(P_mm), PET_mm = sum(PET_mm), .groups = "drop") %>%
  filter(n_clim >= 365)

# 2.1 Check: gridded P vs. on-site rain gauges (Southeast sites) ----------------
gauge_daily <- map_dfr(c("PRF", "TAL", "WHR"), function(ws) {
  fs <- list.files(file.path("aims_data", paste0(ws, "_METS")),
                   pattern = "^METS_.*\\.csv$", full.names = TRUE)
  if (length(fs) == 0) return(NULL)
  map_dfr(fs, ~read_csv(.x, col_types = cols_only(datetime = "c", rainfall = "d"))) %>%
    mutate(date = as_date(ymd_hms(datetime))) %>%
    filter(!is.na(rainfall)) %>%
    group_by(date) %>%
    summarise(gauge_mm = sum(rainfall), .groups = "drop") %>%
    mutate(watershed = ws)
})

p_check <- gauge_daily %>%
  inner_join(climate_daily %>% select(watershed, date, P_mm), by = c("watershed", "date")) %>%
  mutate(WY = water_year(date)) %>%
  group_by(watershed, WY) %>%
  summarise(n_gauge_days = n(), gauge_mm = sum(gauge_mm), gridded_mm = sum(P_mm),
            .groups = "drop") %>%
  filter(n_gauge_days >= 0.9 * 365) %>%      # nearly complete gauge years only
  mutate(ratio_gridded_gauge = round(gridded_mm / gauge_mm, 2))
cat("\ngridMET vs. on-site gauge, water-year totals (days with gauge data):\n")
print(p_check)

# 3.0 Annual runoff and network wetness ----------------------------------------
hydro_wy <- read_csv(merge_file, show_col_types = FALSE) %>%
  filter(watershed %in% site_order) %>%
  mutate(WY = water_year(date)) %>%
  group_by(watershed, WY) %>%
  summarise(days_in_wy = as.numeric(as_date(paste0(first(WY), "-09-30")) -
                                      as_date(paste0(first(WY) - 1, "-10-01"))) + 1,
            n_Q      = sum(!is.na(Q_mm_day)),
            n_ASDN   = sum(!is.na(ASDN_pct)),
            Q_mm     = mean(Q_mm_day, na.rm = TRUE) * days_in_wy,
            ASDN_mean = mean(ASDN_pct, na.rm = TRUE),
            pct_days_flowing = 100 * mean(Q_mm_day > 0, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Q_cov    = n_Q / days_in_wy,
         ASDN_cov = n_ASDN / days_in_wy)

# 4.0 Budyko table ---------------------------------------------------------------
budyko <- hydro_wy %>%
  inner_join(climate_wy, by = c("watershed", "WY")) %>%
  filter(Q_cov >= cov_min) %>%
  mutate(ASDN_mean    = if_else(ASDN_cov >= asdn_cov_min, ASDN_mean, NA_real_),
         completeness = if_else(Q_cov >= cov_complete, "Complete (≥ 85% of days with Q)",
                                "Incomplete (50–85% of days with Q)"),
         aridity      = PET_mm / P_mm,
         runoff_ratio = Q_mm / P_mm,
         evap_index   = 1 - runoff_ratio,
         watershed    = factor(watershed, levels = site_order)) %>%
  left_join(site_info %>% mutate(watershed = factor(watershed, levels = site_order)),
            by = "watershed") %>%
  arrange(watershed, WY)

budyko %>%
  select(watershed, region, WY, Q_cov, ASDN_cov, P_mm, PET_mm, Q_mm, aridity,
         evap_index, ASDN_mean, pct_days_flowing, completeness) %>%
  mutate(across(c(Q_cov, ASDN_cov), ~round(.x, 2)),
         across(c(P_mm, PET_mm, Q_mm, ASDN_mean, pct_days_flowing), ~round(.x, 0)),
         across(c(aridity, evap_index), ~round(.x, 2))) %>%
  write_csv(file.path(tab_dir, "SI_table_budyko_water_years.csv"))
print(select(budyko, watershed, WY, Q_cov, P_mm, PET_mm, Q_mm, aridity, evap_index,
             ASDN_mean), n = Inf)

# 5.0 Figure: Budyko space, one point per watershed-water year ----------------
budyko_curve <- tibble(aridity = 10^seq(log10(0.2), log10(5), length.out = 300)) %>%
  mutate(evap_index = sqrt(aridity * tanh(1 / aridity) * (1 - exp(-aridity))))
limits <- tibble(aridity = c(0, 1, 5), evap_index = c(0, 1, 1))

x_max <- max(1.6, ceiling(max(budyko$aridity) * 10) / 10 + 0.1)
y_min <- min(0.2, floor(min(budyko$evap_index) * 10) / 10 - 0.05)

# Point labels: site code and water year (e.g., "PRF 22")
budyko <- budyko %>% mutate(pt_label = paste(watershed, substr(WY, 3, 4)))

fig_budyko <- ggplot(budyko, aes(x = aridity, y = evap_index)) +
  geom_line(data = limits, color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_line(data = budyko_curve, color = "black", linewidth = 0.7) +
  geom_point(data = filter(budyko, is.na(ASDN_mean)), aes(color = watershed),
             shape = 4, size = 2.5, stroke = 1, show.legend = FALSE) +
  geom_point(data = filter(budyko, !is.na(ASDN_mean)),
             aes(size = ASDN_mean, fill = watershed, alpha = completeness),
             shape = 21, color = "black", stroke = 0.4) +
  geom_text_repel(aes(label = pt_label, color = watershed), size = 2.6, seed = 3,
                  box.padding = 0.6, point.padding = 0.4, min.segment.length = 0.2,
                  bg.color = "white", bg.r = 0.15,
                  segment.color = "grey60", segment.size = 0.3, max.overlaps = Inf,
                  show.legend = FALSE) +
  annotate("text", x = 0.62, y = 0.70, label = "Energy limit", angle = 42, size = 2.8,
           color = "grey40") +
  annotate("text", x = 1.02, y = 1.015, label = "Water limit", hjust = 0,
           size = 2.8, color = "grey40") +
  annotate("text", x = 1.3, y = 0.66, label = "Budyko (1974)", size = 2.8, angle = 20) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_fill_manual(values = col_line, name = "Watershed") +
  scale_alpha_manual(values = c(1, 0.35), name = "Discharge record") +
  scale_size_area(max_size = 9, breaks = c(25, 50, 75, 100), limits = c(0, 100),
                  name = "Mean ASDN over\nwater year (%)") +
  coord_cartesian(xlim = c(0.4, x_max), ylim = c(y_min, 1.05)) +
  guides(fill  = guide_legend(order = 1, override.aes = list(size = 4)),
         size  = guide_legend(order = 2, override.aes = list(fill = "grey70")),
         alpha = guide_legend(order = 3, override.aes = list(size = 4, fill = "grey40"))) +
  labs(x = "Aridity index, PET / P", y = "Evaporative index, (P − Q) / P",
       caption = paste0(
         "One point per watershed and water year (Oct–Sep). P and PET (grass reference ET): ",
         climate_source, " at each outlet. Q: AIMS outlet discharge.\n",
         "Complete: ≥ 85% of days with Q; incomplete: 50–85% (annual runoff uncertain). ",
         "×: too little ASDN data for a wetness value.")) +
  theme_bw() +
  theme(axis.title = element_text(size = 11), panel.grid.minor = element_blank(),
        legend.position = "right", legend.title = element_text(size = 9),
        legend.text = element_text(size = 8), plot.caption = element_text(size = 8, hjust = 0))
fig_budyko

ggsave(file.path(fig_dir, "SI_budyko_ASDN_water_years.png"), fig_budyko,
       width = 8, height = 6.2, dpi = 300)

# 6.0 Figure: Budyko space, one point per watershed ----------------------------
# Point = mean of the annual indices over the water years shown in section 5
# (all years with >= 50% Q coverage). Bars = range (min-max) across those years.
# Size = mean of the annual mean ASDN. With 2-3 years per site the range is a
# better description of interannual variability than an SD.
budyko_site <- budyko %>%
  group_by(watershed, region) %>%
  summarise(n_wy          = n(),
            n_complete    = sum(Q_cov >= cov_complete),
            wy_list       = paste(WY, collapse = ", "),
            aridity_mean  = mean(aridity),
            aridity_min   = min(aridity),
            aridity_max   = max(aridity),
            evap_mean     = mean(evap_index),
            evap_min      = min(evap_index),
            evap_max      = max(evap_index),
            ASDN_mean     = mean(ASDN_mean, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(record = if_else(n_complete == n_wy, "All years complete",
                          "Includes incomplete years"),
         site_label = sprintf("%s (n = %d)", watershed, n_wy))

budyko_site %>%
  mutate(across(where(is.double), ~round(.x, 2))) %>%
  write_csv(file.path(tab_dir, "SI_table_budyko_site_means.csv"))
print(budyko_site, width = Inf)

fig_budyko_site <- ggplot(budyko_site, aes(x = aridity_mean, y = evap_mean)) +
  geom_line(data = limits, aes(x = aridity, y = evap_index),
            color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_line(data = budyko_curve, aes(x = aridity, y = evap_index),
            color = "black", linewidth = 0.7) +
  geom_errorbar(aes(ymin = evap_min, ymax = evap_max, color = watershed),
                width = 0, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = aridity_min, xmax = aridity_max, color = watershed),
                 height = 0, linewidth = 0.7) +
  geom_point(aes(size = ASDN_mean, fill = watershed, alpha = record),
             shape = 21, color = "black", stroke = 0.5) +
  geom_text_repel(aes(label = site_label, color = watershed), size = 3.2,
                  fontface = "bold", seed = 4, box.padding = 1, point.padding = 0.6,
                  min.segment.length = 0.3, segment.color = "grey60", segment.size = 0.3,
                  bg.color = "white", bg.r = 0.15, show.legend = FALSE) +
  annotate("text", x = 0.62, y = 0.70, label = "Energy limit", angle = 42, size = 2.8,
           color = "grey40") +
  annotate("text", x = 1.02, y = 1.015, label = "Water limit", hjust = 0,
           size = 2.8, color = "grey40") +
  annotate("text", x = 2.45, y = 0.965, label = "Budyko (1974)", size = 2.8, angle = 6, hjust = 1) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_fill_manual(values = col_line, guide = "none") +
  scale_alpha_manual(values = c("All years complete" = 1, "Includes incomplete years" = 0.4),
                     name = "Discharge record") +
  scale_size_area(max_size = 10, breaks = c(25, 50, 75, 100), limits = c(0, 100),
                  name = "Mean ASDN (%)") +
  coord_cartesian(xlim = c(0.4, x_max), ylim = c(y_min, 1.05)) +
  guides(size  = guide_legend(order = 1, override.aes = list(fill = "grey70")),
         alpha = guide_legend(order = 2, override.aes = list(size = 4, fill = "grey40"))) +
  labs(x = "Aridity index, PET / P", y = "Evaporative index, (P − Q) / P",
       caption = paste0(
         "Points: mean of annual values over the water years with ≥ 50% of days with Q ",
         "(n = number of years). Bars: range across years.\n",
         "Point size: mean network extent (ASDN, % of maximum active length). P and PET (grass ",
         "reference ET): ", climate_source, ". Q: AIMS outlet discharge.")) +
  theme_bw() +
  theme(axis.title = element_text(size = 11), panel.grid.minor = element_blank(),
        legend.position = "right", legend.title = element_text(size = 9),
        legend.text = element_text(size = 8), plot.caption = element_text(size = 8, hjust = 0))
fig_budyko_site

ggsave(file.path(fig_dir, "SI_budyko_ASDN.png"), fig_budyko_site,
       width = 8, height = 6.2, dpi = 300)
