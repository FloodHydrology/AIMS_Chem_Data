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

# 4.0 Quality assurance checks (region-agnostic) -------------------------------
# Portable across regions: no hardcoded sites/dates. Collects every issue into
# qa_log and writes it per region. Flags issues; applies only rules you've OK'd.
region_tag <- unique(na.omit(chem_wide$watershed)) |> paste(collapse = "_")
qa_log <- list()

# 4.1 Watershed mapping integrity — THE portability check.
# recode() silently yields NA for unmapped site prefixes. In a new region that
# means every row. Fail loud here instead of exporting NA-watershed data.
unmapped <- chem_wide |> filter(is.na(watershed))
qa_log$unmapped_watershed <- nrow(unmapped)
if (nrow(unmapped) > 0) {
  warning("QA 4.1: ", nrow(unmapped), " rows have NA watershed — site prefix not in recode(). ",
          "Update the watershed map for this region before trusting output.")
  print(unmapped |> distinct(siteId))
}

# 4.2 Duplicate date+site rows (source reanalyses; fan out downstream joins)
dupes <- chem_wide |> group_by(date, siteId, watershed) |> filter(n() > 1) |> ungroup()
qa_log$duplicate_rows <- nrow(dupes)
cat("\n[4.2] Duplicate date+site rows:", nrow(dupes), "across",
    n_distinct(dupes[c("date","siteId")]), "combos\n")
if (nrow(dupes)) print(arrange(dupes, siteId, date) |> select(date, siteId, watershed))

# 4.3 Date range — catch a mis-parsed date (parser handles 3 encodings)
qa_log$na_date <- sum(is.na(chem_wide$date))
qa_log$out_of_window <- chem_wide |>
  filter(date < as.Date("2021-07-01") | date > as.Date("2024-12-31")) |> nrow()
cat("\n[4.3] NA dates:", qa_log$na_date, "| out-of-window:", qa_log$out_of_window, "\n")

# 4.4 Non-positive analyte values — impossible concentrations; log-plot drops
# them silently. RULE (applied): set <=0 to NA, logged so nothing vanishes unseen.
analyte_cols <- c(anio_cols, caio_cols, nutr_cols, docs_cols)
nonpos <- chem_wide |>
  pivot_longer(all_of(analyte_cols), names_to="analyte", values_to="value") |>
  filter(!is.na(value) & value <= 0)
qa_log$nonpositive <- nrow(nonpos)
cat("\n[4.4] Non-positive values set to NA:", nrow(nonpos), "\n")
if (nrow(nonpos)) print(count(nonpos, analyte, name="n"))
chem_wide <- chem_wide |> mutate(across(all_of(analyte_cols), ~ ifelse(.x <= 0, NA, .x)))

# 4.5 Empty rows — joined-in date/site with no analyte data at all
qa_log$empty_rows <- chem_wide |> filter(if_all(all_of(analyte_cols), is.na)) |> nrow()
cat("\n[4.5] All-NA analyte rows:", qa_log$empty_rows, "\n")

# 4.6 No approach flag — ENVI roster / coalesce sanity
qa_log$no_approach <- chem_wide |> filter(appr1==0&appr2==0&appr3==0&appr4==0) |> nrow()
cat("\n[4.6] Rows with no approach flag:", qa_log$no_approach, "\n")

# 4.7 Completeness per analyte (expected sparsity varies by region)
cat("\n[4.7] NA per analyte (of", nrow(chem_wide), "rows):\n")
print(chem_wide |> summarise(across(all_of(analyte_cols), ~ sum(is.na(.x)))) |>
        pivot_longer(everything(), names_to="analyte", values_to="n_NA") |>
        arrange(desc(n_NA)) |> as.data.frame())

# Write per-region QA record so each region leaves an auditable trail
qa_summary <- tibble(region = region_tag, check = names(qa_log),
                     value = unlist(qa_log))
write_csv(qa_summary, paste0("output/qa_log_", region_tag, ".csv"))
cat("\nQA log written: output/qa_log_", region_tag, ".csv\n", sep="")
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