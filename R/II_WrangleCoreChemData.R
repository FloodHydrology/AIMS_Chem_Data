#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Wrangle core chem data
# Date: 8/5/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Wrangle core chem datasets at outlet + LTM sites of AIMS watersheds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.0 Setup workspace ----------------------------------------------------------
#Clear workspace
remove(list=ls())

#install libraries of interest
library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)

# 2.0 Wrangle data -------------------------------------------------------------
#Download core datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
files<-list.files(
  "aims_data", 
  pattern = "(NUTR|ANIO|CAIO|DOCS).*\\.xlsx$",
  recursive = TRUE, full.names = TRUE)

#download data
all_chem_data <- setNames(lapply(files, read_excel, sheet = "Final Data"), basename(files))

#Harmonize site column name (ANIO files use "Site", others use "siteId")
all_chem_data <- map(all_chem_data, ~ rename(.x, siteId = any_of(c("Site", "siteId"))))

#Identify outlet + LTM sites from ENVI metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ENVI uses appr1 (outlet) / appr2 (long-term monitoring) to tag SITE ROLE.
# Sites with neither flag are STIC/synoptic-only and are excluded.
envi_files <- list.files("aims_data", pattern = "ENVI.*\\.xlsx$",
                         recursive = TRUE, full.names = TRUE)

ltm_sites <- envi_files |>
  map(~ read_excel(.x, sheet = "Final Data") |>
        rename(siteId = any_of(c("Site", "siteId"))) |>
        mutate(across(any_of(c("appr1","appr2")), as.numeric)) |>
        filter(appr1 == 1 | appr2 == 1) |>
        pull(siteId)) |>
  unlist() |> unique()

cat("Outlet + LTM sites:", length(ltm_sites), "\n"); print(sort(ltm_sites))

#Filter chem data to outlet + LTM sites (all rows for those sites, any approach)
outlet_chem_data <- map(all_chem_data, ~ filter(.x, siteId %in% ltm_sites))

#fix inconsistent date formats~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Harmonize the date column name (ANIO uses "Date", others "date")
outlet_chem_data <- map(outlet_chem_data, ~ rename(.x, date = any_of(c("Date", "date"))))

# Robust parser: handles YYYYMMDD ints, Excel serials, and real datetimes
#    (branches on magnitude — YYYYMMDD ~2e7, Excel serial ~4e4)
parse_mixed_date <- function(x) {
  if (inherits(x, "POSIXct") || inherits(x, "Date")) return(as.Date(x))
  x <- as.numeric(x)
  out <- as.Date(ifelse(
    x > 1e7,
    as.numeric(as.Date(as.character(x), format = "%Y%m%d")),
    as.numeric(as.Date(x, origin = "1899-12-30"))
  ), origin = "1970-01-01")
  out
}

# Apply across all 12
outlet_chem_data <- map(outlet_chem_data, ~ mutate(.x, date = parse_mixed_date(date)))

# verify parser actually worked
# every date should now fall inside the study window; NAs = parse failures
map(outlet_chem_data, ~ range(.x$date, na.rm = TRUE))
map_int(outlet_chem_data, ~ sum(is.na(.x$date)))   # want all zeros

#Add watershed column from the site-code prefix (TL->TAL, PR->PRF, WH->WHR)~~~~
outlet_chem_data <- map(outlet_chem_data, ~ .x |>
                          mutate(watershed = recode(substr(siteId, 1, 2),
                                                    TL = "TAL", PR = "PRF", WH = "WHR")))

#Create a single wide format csv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine the watersheds within each data type, keeping date/siteId + analytes
# NOTE: coerce analytes to numeric PER FILE before binding, so character/double
#       clashes across files (e.g. Sulfate_avg) don't break bind_rows()
# NOTE: approach flags (appr1-4) carried through too; ANIO lacks appr4 so any_of()
#       tolerates it, and fix_appr() below fills it as 0
appr_cols <- c("appr1","appr2","appr3","appr4")

bind_type <- function(keys, analytes) {
  keys |>
    map(~ .x |>
          select(date, siteId, watershed, any_of(appr_cols), any_of(analytes)) |>
          mutate(across(any_of(analytes), as.numeric))) |>  # coerce first
    bind_rows()                                              # then stack
}
anio_cols <- c("Flouride_avg","Chloride_avg","Nitrite_avg","Bromide_avg",
               "Nitrate_avg","Phosphate_avg","Sulfate_avg")
caio_cols <- c("Na_avg","Ca_avg","B_avg","Mg_avg","Si_avg","K_avg","Sr_avg")
nutr_cols <- c("SRPugL","NH4ugL","NO3NO2ugL","NO3ugL")
docs_cols <- c("DOCmgL")
anio <- bind_type(outlet_chem_data[grep("^ANIO", names(outlet_chem_data))], anio_cols)
caio <- bind_type(outlet_chem_data[grep("^CAIO", names(outlet_chem_data))], caio_cols)
nutr <- bind_type(outlet_chem_data[grep("^NUTR", names(outlet_chem_data))], nutr_cols)
docs <- bind_type(outlet_chem_data[grep("^DOCS", names(outlet_chem_data))], docs_cols)

# normalize approach flags: ANIO has no appr4 -> add it; NA -> 0; make integer
fix_appr <- function(df) {
  if (!"appr4" %in% names(df)) df$appr4 <- 0
  df |> mutate(across(all_of(appr_cols),
                      ~ as.integer(replace_na(as.numeric(.x), 0))))
}
anio <- fix_appr(anio); caio <- fix_appr(caio)
nutr <- fix_appr(nutr); docs <- fix_appr(docs)

# Pull approach flags out to their own per-sample table, coalesced across the four
# data types (max = if any type tagged this date/site with an approach, keep it).
# This avoids appr1.x/appr1.y collisions and keeps approach flags off the join key,
# so differing flags across types can't fragment the table into duplicate rows.
appr_tbl <- bind_rows(anio, caio, nutr, docs) |>
  group_by(date, siteId, watershed) |>
  summarise(across(all_of(appr_cols), ~ max(.x, na.rm = TRUE)), .groups = "drop")

# drop the per-type approach flags before joining analytes (re-added once, below)
anio <- select(anio, -all_of(appr_cols))
caio <- select(caio, -all_of(appr_cols))
nutr <- select(nutr, -all_of(appr_cols))
docs <- select(docs, -all_of(appr_cols))

# join all four analyte sets, then attach the single coalesced approach table
chem_wide <- list(anio, caio, nutr, docs) |>
  reduce(full_join, by = c("date", "siteId", "watershed")) |>
  left_join(appr_tbl, by = c("date", "siteId", "watershed")) |>
  relocate(all_of(appr_cols), .after = watershed) |>
  arrange(siteId, date)

# 4.0 Flag duplicates ----------------------------------------------------------
# PLACEHOLDER — next step. Find date+siteId combos with >1 row (true source
# duplicates / reanalyses), e.g. TLM01 2022-06-09 (ANIO) and 2023-01-30 (DOCS).
# Expect more now that LTM sites are included.

# 5.0 Export -------------------------------------------------------------------
#Plots plots plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# long form for faceting (exclude id + approach cols from the pivot)
plot_df <- chem_wide |>
  pivot_longer(
    cols = -c(date, siteId, watershed, all_of(appr_cols)),
    names_to = "analyte", values_to = "value"
  ) |>
  filter(!is.na(value), value > 0)              # drop NA + non-positive for log scale

ggplot(plot_df, aes(x = date, y = value, color = watershed)) +
  geom_line(lwd = 0.5, linetype = "dashed") +
  geom_point(size = 1.2) +
  scale_y_log10() +
  scale_color_manual(values = c(TAL = "steelblue", PRF = "#B4664E", WHR = "#5B8C5A"),
                     name = "Watershed") +
  facet_wrap(~ analyte, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 8),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "top"
  ) +
  xlab(NULL) +
  ylab("Concentration")

#Export csv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write_csv(chem_wide, "output/outlet_chem_data.csv")