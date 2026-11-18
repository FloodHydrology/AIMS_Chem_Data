#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Hydrologic analysis: outlet discharge vs. network connectivity (ASDN)
# Date: 11/18/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Describe how active stream drainage network extent (ASDN) scales with
#          outlet streamflow (Q) across the AIMS watersheds: elasticity (slope of
#          the log-log power law) and dispersion (scatter around it). Hydrologic
#          basis for the concentration-connectivity analysis.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sections and outputs (figures -> figures/, tables -> output/)
#  1.0 Setup and data
#  2.0 Figure    Time series of Q and ASDN, all sites
#  3.0 Figure    Duration curves of Q and ASDN with sampling days
#      Table     Where sampling days fall on the Q and ASDN duration curves
#  4.0 Figure    log-log ASDN vs. Q, power-law fit per site
#      Table     Power-law fits: elasticity, R2, days used/excluded
#  5.0 Figure    Residuals of the power-law fits vs. Q ("dispersion")
#  6.0 Figure    ASDN vs. Q colored by season
#      Figure    Residuals vs. Q colored by season
#      Figure    Residuals vs. Q colored by day of year
#      Figure    Residuals vs. antecedent wetness
#      Table     Dispersion and the R2 added by limb, season, antecedent wetness
#  7.0 Main fig 1 WHR: time series, elasticity, dispersion, STIC maps on two days
#                (uses R/spatial_functions.R and aims_data/WHR_ENVI, WHR_STIC)
#  8.0 Main fig 2 Cross-site: ASDN vs. Q with power laws, elasticity vs. dispersion
#
# Inputs
#  data/AIMS_Q_connectivity_merge/output/streamDAG_Q_data.csv
#     Merged daily streamDAG connectivity metrics + mean daily outlet discharge,
#     one row per watershed-day (E. Seybold; see that folder's README). Read as-is.
#  data/allSites_chem_Q_data.csv   Outlet chemistry sampling days (is_outlet)
#  data/site_colors.xlsx           Site colors: hex2 = lines/points, hex = fills
#
# Key decisions
#  - Sites: KNZ, YMR (Great Plains), GBJ, JDR (Mountain West), PRF, TAL, WHR
#    (Southeast). DRC excluded (no discharge yet).
#  - ASDN = ASDN_pct: active network length as % of the site's maximum active
#    length over the full streamDAG record.
#  - Q = Q_mm_day (area-normalized). GBJ area is unconfirmed (16.2 vs 14.83 km2);
#    this shifts GBJ along the Q axis but does not change its elasticity.
#  - Q floor = 0.001 mm/day. Days with 0 < Q < floor (down to 1e-18 mm/day at PRF,
#    likely rating-curve artifacts near zero stage) are treated as zero flow in the
#    fits and drawn at the floor in plots. Without the floor PRF elasticity drops
#    from 0.22 to 0.09.
#  - Exceedance = % of days with a value >= that day's value. Tied days share the
#    exceedance of their block's right edge (ASDN has long tied runs at 100%).
#  - Power law fit by OLS on log10 values per site; days with Q < floor or
#    ASDN = 0 excluded. Dispersion = SD of the log10 residuals.
#
# Open data questions (for E. Seybold)
#  - PRF: ~16% ASDN band at high Q, mostly Aug-Nov 2021 (sensor deployment?).
#  - KNZ, YMR: low ASDN in winter at moderate Q. Are frozen STICs read as dry?
#  - JDR (265 days) and YMR (41 days) have streamDAG rows with ASDN_pct = NA.
#  - How is each node's daily wet/dry state defined, and are gaps filled?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1.0 Setup workspace ----------------------------------------------------------
#Clear workspace
remove(list = ls())

#Load libraries of interest
library(tidyverse)
library(lubridate)
library(readxl)
library(patchwork)

#Define file paths
merge_file  <- "data/AIMS_Q_connectivity_merge/output/streamDAG_Q_data.csv"
chem_file   <- "data/allSites_chem_Q_data.csv"
colors_file <- "data/site_colors.xlsx"
fig_dir     <- "figures"
tab_dir     <- "output"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)

#Analysis settings
Q_floor       <- 1e-3   # mm/day; lower Q limit for fits and log axes
limb_tol      <- 0.01   # fractional day-to-day change in Q that counts as rising/falling
ante_days     <- 30     # window (days) for antecedent wetness
ante_min_days <- 20     # minimum days with Q in that window

# 1.1 Sites, order, and colors -------------------------------------------------
site_info <- tribble(
  ~watershed, ~region,
  "KNZ",      "Great Plains",
  "YMR",      "Great Plains",
  "GBJ",      "Mountain West",
  "JDR",      "Mountain West",
  "PRF",      "Southeast",
  "TAL",      "Southeast",
  "WHR",      "Southeast"
)
site_order <- site_info$watershed
site_title <- setNames(paste0(site_info$watershed, " (", site_info$region, ")"),
                       site_info$watershed)
site_lab   <- as_labeller(site_title)

site_colors <- read_excel(colors_file) %>% rename(watershed = site)
stopifnot(all(site_order %in% site_colors$watershed))

col_line <- deframe(select(site_colors, watershed, hex2))[site_order]   # lines, points
col_fill <- deframe(select(site_colors, watershed, hex))[site_order]    # fills

# 1.2 Daily Q and ASDN ---------------------------------------------------------
hydro <- read_csv(merge_file, show_col_types = FALSE) %>%
  filter(watershed %in% site_order) %>%
  mutate(watershed = factor(watershed, levels = site_order)) %>%
  select(watershed, siteId, date, ASDN_pct, Q_Ls, Q_mm_day, Q_flag, Q_partial_flag) %>%
  arrange(watershed, date)

# Check: record length, days with each variable, and zeros
hydro %>%
  group_by(watershed) %>%
  summarise(start       = min(date),
            end         = max(date),
            n_days      = n(),
            n_both      = sum(!is.na(Q_mm_day) & !is.na(ASDN_pct)),
            n_Q_zero    = sum(Q_mm_day == 0, na.rm = TRUE),
            n_ASDN_zero = sum(ASDN_pct == 0, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# 1.3 Chemistry sampling days --------------------------------------------------
samples <- read_csv(chem_file, show_col_types = FALSE) %>%
  filter(is_outlet, watershed %in% site_order) %>%
  transmute(watershed = factor(watershed, levels = site_order),
            date      = as_date(date)) %>%
  distinct()

# 1.4 Shared plot settings -----------------------------------------------------
theme_hydro <- theme_bw() +
  theme(axis.title       = element_text(size = 12),
        axis.text        = element_text(size = 9),
        strip.text       = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank())

log_labels <- scales::label_number(drop0trailing = TRUE, big.mark = ",")

# Panel title in the site color (used by the patchwork figures)
site_title_theme <- function(ws) {
  theme(plot.title  = element_text(color = col_line[ws], face = "bold", size = 11,
                                   margin = margin(b = 2)),
        plot.margin = margin(2, 6, 2, 4))
}

# 2.0 Figure: time series of Q and ASDN ----------------------------------------
# One panel per site on a shared date axis (style after Peterson et al.):
#   ASDN on the left axis (site color), Q on a log right axis (gray),
#   ticks along the bottom on chemistry sampling days.
# Each site needs its own Q-axis transform, so panels are built separately.

ts_data <- hydro %>%
  group_by(watershed) %>%
  complete(date = seq(min(date), max(date), by = "day")) %>%    # gaps break lines
  ungroup() %>%
  mutate(Q_plot = case_when(Q_mm_day <= 0 ~ NA_real_,             # zero flow -> gap
                            TRUE          ~ pmax(Q_mm_day, Q_floor)))
date_lim <- range(ts_data$date)

plot_ts_site <- function(ws, letter, show_x = FALSE, show_ylab = FALSE) {
  d <- filter(ts_data, watershed == ws)
  s <- filter(samples, watershed == ws)

  # Map log10(Q) onto the 0-100 ASDN scale using whole decades of this site's Q
  q_lo <- floor(log10(min(d$Q_plot, na.rm = TRUE)))
  q_hi <- ceiling(log10(max(d$Q_plot, na.rm = TRUE)))
  to_y <- function(q) 100 * (log10(q) - q_lo) / (q_hi - q_lo)
  to_q <- function(y) 10^(q_lo + y * (q_hi - q_lo) / 100)

  ggplot(d, aes(x = date)) +
    geom_line(aes(y = to_y(Q_plot)), color = "grey25", linewidth = 0.35) +
    geom_line(aes(y = ASDN_pct), color = col_line[ws], linewidth = 0.7) +
    geom_rug(data = s, aes(x = date), sides = "b", inherit.aes = FALSE,
             color = "grey20", linewidth = 0.5, length = unit(0.06, "npc")) +
    annotate("text", x = date_lim[2], y = 100, label = letter,
             hjust = 1.2, vjust = 1.1, fontface = "bold", size = 5) +
    scale_y_continuous(
      name     = if (show_ylab) "ASDN (% of maximum active length)" else NULL,
      limits   = c(0, 100), breaks = seq(0, 100, 25),
      sec.axis = sec_axis(transform = to_q, breaks = 10^(q_lo:q_hi), labels = log_labels,
                          name = if (show_ylab) "Discharge (mm/day)" else NULL)) +
    scale_x_date(limits = date_lim, date_breaks = "6 months", date_labels = "%Y-%m",
                 expand = expansion(mult = 0.01)) +
    labs(title = site_title[[ws]], x = if (show_x) "Date" else NULL) +
    theme_hydro +
    site_title_theme(ws) +
    theme(axis.text.x  = if (show_x) element_text(size = 9) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.text.y  = element_text(size = 8))
}

fig_timeseries <- map2(site_order, LETTERS[seq_along(site_order)], function(ws, l) {
  plot_ts_site(ws, l,
               show_x    = ws == last(site_order),
               show_ylab = ws == site_order[ceiling(length(site_order) / 2)])
}) %>%
  wrap_plots(ncol = 1) +
  plot_annotation(caption = paste0(
    "Colored line: ASDN (left axis). Gray line: daily mean discharge (right axis, log scale; ",
    "gaps = zero flow or no data; values < ", Q_floor, " mm/day drawn at ", Q_floor, ").\n",
    "Ticks: outlet chemistry sampling days."))
fig_timeseries

ggsave(file.path(fig_dir, "SI_timeseries_Q_ASDN.png"), fig_timeseries,
       width = 10, height = 13, dpi = 300)

# 3.0 Figure and table: duration curves with sampling days ------------------
# For each site and variable:
#   exceed_pct = % of days with a value >= that day's value (daily record)
#   curve      = step function from x = 0 (site maximum) to x = 100
#   vertical lines = sampling days; a number = sampling days sharing one line
#   shading    = block of days with a value of zero

duration_data <- function(var) {
  hydro %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(watershed) %>%
    mutate(value      = .data[[var]],
           exceed_pct = 100 * rank(-value, ties.method = "max") / n()) %>%
    ungroup() %>%
    select(watershed, date, value, exceed_pct)
}

# Step-curve vertices: one row per distinct value, plus x = 0 at the maximum
duration_curve <- function(d) {
  d <- d %>% distinct(value, exceed_pct) %>% arrange(exceed_pct)
  bind_rows(tibble(value = d$value[1], exceed_pct = 0), d)
}

dur <- list(
  Q    = list(var = "Q_mm_day", log_y = TRUE,  ylab = "Daily mean discharge (mm/day)"),
  ASDN = list(var = "ASDN_pct", log_y = FALSE, ylab = "ASDN (% of maximum active length)")
)

for (k in names(dur)) {
  dc <- duration_data(dur[[k]]$var)
  dur[[k]]$data    <- dc
  dur[[k]]$samples <- samples %>% inner_join(dc, by = c("watershed", "date"))
  dur[[k]]$zeros   <- dc %>%
    group_by(watershed) %>%
    summarise(xmin = 100 * sum(value > 0) / n(),
              pct_zero = round(100 * mean(value == 0)), .groups = "drop") %>%
    filter(pct_zero > 0)
}

duration_panel <- function(k, ws, letter, title = NULL, subtitle = NULL) {
  log_y <- dur[[k]]$log_y
  d <- dur[[k]]$data %>%
    filter(watershed == ws, if (log_y) value > 0 else TRUE) %>%
    duration_curve()
  z <- filter(dur[[k]]$zeros, watershed == ws)
  s <- filter(dur[[k]]$samples, watershed == ws)
  s_stack <- s %>% count(exceed_pct) %>% filter(n > 1)

  # Panel bottom: on a log axis use 0 (-> -Inf after transform), not -Inf (-> NaN)
  y_bottom <- if (log_y) 0 else -Inf
  y_label  <- if (log_y) 0 else 5

  p <- ggplot() +
    geom_rect(data = z, aes(xmin = xmin, xmax = 100, ymin = y_bottom, ymax = Inf),
              fill = col_fill[ws], alpha = 0.35) +
    geom_text(data = filter(z, pct_zero >= 5),
              aes(x = (xmin + 100) / 2, y = y_label, label = paste0("zero\n", pct_zero, "%")),
              vjust = -0.4, size = 2.6, color = "grey20") +
    geom_vline(data = s, aes(xintercept = exceed_pct),
               color = "grey30", alpha = 0.6, linewidth = 0.3) +
    geom_text(data = s_stack,
              aes(x = exceed_pct, y = Inf, label = n,
                  hjust = if_else(exceed_pct > 90, 1.3, -0.2),
                  vjust = if_else(exceed_pct > 90, 3.2, 1.3)),
              size = 2.4, color = "grey20") +
    geom_step(data = d, aes(x = exceed_pct, y = value), direction = "vh",
              color = col_line[ws], linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = letter,
             hjust = 1.5, vjust = 1.4, fontface = "bold", size = 4) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = 0.01)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_hydro +
    site_title_theme(ws) +
    theme(axis.text = element_text(size = 8),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5))

  if (log_y) {
    p + scale_y_log10(limits = c(Q_floor, NA), oob = scales::squish, labels = log_labels)
  } else {
    p + scale_y_continuous(limits = c(0, 100))
  }
}

# Q in the left column, ASDN in the right, one row per site
fig_duration <- map(seq_along(site_order), function(i) {
  ws <- site_order[i]
  list(duration_panel("Q", ws, LETTERS[2 * i - 1], title = site_title[[ws]],
                      subtitle = if (i == 1) dur$Q$ylab),
       duration_panel("ASDN", ws, LETTERS[2 * i], title = " ",
                      subtitle = if (i == 1) dur$ASDN$ylab))
}) %>%
  flatten() %>%
  wrap_plots(ncol = 2) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(caption = paste0(
    "Vertical lines: outlet chemistry sampling days (numbers: sampling days sharing one line). ",
    "Shading: days with a value of zero.\nExceedance from the daily record. ",
    "Discharge values < ", Q_floor, " mm/day drawn at ", Q_floor, ".")) &
  labs(x = "Percent of days equaled or exceeded")
fig_duration

ggsave(file.path(fig_dir, "SI_duration_curves_samples.png"), fig_duration,
       width = 10, height = 16, dpi = 300)

# Table: where sampling days fall on each duration curve
table_sampling <- map_dfr(names(dur), function(k) {
  dur[[k]]$samples %>%
    group_by(watershed) %>%
    summarise(n_samples       = n(),
              exceed_median   = median(exceed_pct),
              exceed_q25      = quantile(exceed_pct, 0.25),
              exceed_q75      = quantile(exceed_pct, 0.75),
              exceed_min      = min(exceed_pct),
              exceed_max      = max(exceed_pct),
              n_samples_zero  = sum(value == 0),
              .groups = "drop") %>%
    left_join(dur[[k]]$data %>%
                group_by(watershed) %>%
                summarise(pct_days_zero = 100 * mean(value == 0), .groups = "drop"),
              by = "watershed") %>%
    mutate(variable = k, .after = watershed)
}) %>%
  arrange(watershed, variable) %>%
  mutate(across(where(is.double), ~round(.x, 1)))
print(table_sampling, n = Inf)

write_csv(table_sampling, file.path(tab_dir, "SI_table_sampling_coverage.csv"))

# 4.0 Figure and table: power-law fits --------------------------------------
# ASDN = a * Q^b, fit by OLS in log10 space for each site:
#   log10(ASDN_pct) = log10(a) + b * log10(Q_mm_day)
#   b = elasticity (% change in ASDN per 1% change in Q)
# Days used: Q >= Q_floor and ASDN > 0. ASDN = 100% days are kept and marked.

fit_data <- hydro %>%
  filter(!is.na(Q_mm_day), !is.na(ASDN_pct)) %>%
  mutate(use_in_fit = Q_mm_day >= Q_floor & ASDN_pct > 0,
         lQ         = log10(Q_mm_day),
         lASDN      = log10(ASDN_pct),
         saturated  = ASDN_pct == 100)

fit_one <- function(d) {
  m  <- lm(lASDN ~ lQ, data = d)
  ci <- confint(m)
  tibble(n           = nobs(m),
         elasticity  = coef(m)[["lQ"]],
         b_ci_low    = ci["lQ", 1],
         b_ci_high   = ci["lQ", 2],
         a           = 10^coef(m)[["(Intercept)"]],   # ASDN (%) at Q = 1 mm/day
         r2          = summary(m)$r.squared,
         n_saturated = sum(d$saturated),
         model       = list(m))
}

pl_fits <- fit_data %>%
  filter(use_in_fit) %>%
  group_by(watershed) %>%
  group_modify(~fit_one(.x)) %>%
  ungroup()

# Residuals (log10 observed - log10 predicted) on fit days
fit_data <- fit_data %>%
  left_join(pl_fits %>% select(watershed, model), by = "watershed") %>%
  mutate(resid = if_else(use_in_fit,
                         lASDN - map2_dbl(model, lQ, ~predict(.x, newdata = tibble(lQ = .y))),
                         NA_real_)) %>%
  select(-model)

# Table: power-law fits
table_fits <- pl_fits %>%
  select(-model) %>%
  left_join(fit_data %>%
              group_by(watershed) %>%
              summarise(n_days_Q_and_ASDN = n(),
                        n_excl_Q_zero     = sum(Q_mm_day == 0),
                        n_excl_Q_lt_floor = sum(Q_mm_day > 0 & Q_mm_day < Q_floor),
                        n_excl_ASDN_zero  = sum(ASDN_pct == 0 & Q_mm_day >= Q_floor),
                        .groups = "drop"),
            by = "watershed") %>%
  left_join(site_info %>% mutate(watershed = factor(watershed, levels = site_order)),
            by = "watershed") %>%
  relocate(region, .after = watershed) %>%
  mutate(across(c(elasticity, b_ci_low, b_ci_high, r2), ~round(.x, 3)),
         a = round(a, 1))
print(table_fits, width = Inf)

write_csv(table_fits, file.path(tab_dir, "SI_table_powerlaw_fits.csv"))

# Fitted lines over each site's fitted Q range
pl_lines <- fit_data %>%
  filter(use_in_fit) %>%
  group_by(watershed) %>%
  reframe(lQ = seq(min(lQ), max(lQ), length.out = 50)) %>%
  left_join(pl_fits %>% select(watershed, elasticity, a), by = "watershed") %>%
  mutate(Q = 10^lQ, ASDN_hat = a * Q^elasticity)

fig_powerlaw <- fit_data %>%
  filter(use_in_fit) %>%
  ggplot(aes(x = Q_mm_day, y = ASDN_pct, color = watershed)) +
  geom_point(aes(shape = saturated), size = 1, alpha = 0.4) +
  geom_line(data = pl_lines, aes(x = Q, y = ASDN_hat), color = "black", linewidth = 0.8) +
  geom_text(data = pl_fits,
            aes(x = 0, y = 0, label = sprintf("b = %.2f\nR² = %.2f\nn = %d", elasticity, r2, n)),
            hjust = -0.1, vjust = -0.15, size = 3, color = "black", inherit.aes = FALSE) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c(`FALSE` = "ASDN < 100%", `TRUE` = "ASDN = 100%"), name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = "Daily mean discharge (mm/day)", y = "ASDN (% of maximum active length)",
       caption = paste0("Power law ASDN = a Q^b fit in log10 space per site (black line); ",
                        "days with Q < ", Q_floor, " mm/day or ASDN = 0 excluded.")) +
  theme_hydro +
  theme(legend.position = "bottom")
fig_powerlaw

ggsave(file.path(fig_dir, "SI_ASDN_Q_powerlaw.png"), fig_powerlaw,
       width = 10, height = 9, dpi = 300)

# 5.0 Figure: residuals vs. Q (dispersion) -------------------------------------
# resid > 0: network more extensive than Q predicts; < 0: less extensive.
# Dispersion = residual SD (log10 units); 10^SD = typical multiplicative scatter.

resid_data <- fit_data %>% filter(use_in_fit)

dispersion <- resid_data %>%
  group_by(watershed) %>%
  summarise(resid_sd  = sd(resid),
            factor_sd = 10^resid_sd,
            .groups = "drop")

fig_residuals <- ggplot(resid_data, aes(x = Q_mm_day, y = resid)) +
  geom_rect(data = dispersion, aes(xmin = 0, xmax = Inf, ymin = -resid_sd, ymax = resid_sd),
            fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = watershed, shape = saturated), size = 1, alpha = 0.4) +
  geom_smooth(method = "loess", formula = y ~ x, span = 0.75, se = FALSE,
              color = "black", linewidth = 0.7) +
  geom_text(data = dispersion,
            aes(x = 0, y = Inf, label = sprintf("SD = %.3f\n(×%.2f)", resid_sd, factor_sd)),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  scale_x_log10(labels = log_labels) +
  scale_color_manual(values = col_line, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c(`FALSE` = "ASDN < 100%", `TRUE` = "ASDN = 100%"), name = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = "Daily mean discharge (mm/day)",
       y = "Residual, log10(ASDN observed / ASDN predicted)",
       caption = paste("Residuals from the per-site power-law fits (power-law figure). Gray band: ±1 SD.",
                       "Black line: loess smooth. × = equivalent multiplicative scatter.")) +
  theme_hydro +
  theme(legend.position = "bottom")
fig_residuals

ggsave(file.path(fig_dir, "SI_residuals_vs_Q.png"), fig_residuals,
       width = 10, height = 9, dpi = 300)

# 6.0 Drivers of dispersion: limb, season, antecedent wetness ------------------
# Limb:   day-to-day change in daily Q; rising/falling if > limb_tol, else steady
#         (also steady if the previous day is missing or zero)
# Season: meteorological (DJF, MAM, JJA, SON)
# AWI:    antecedent wetness index = log10(Q_ante / Q_today), Q_ante = mean daily Q
#         over the previous ante_days days (today excluded; >= ante_min_days with
#         data; zeros count). Both terms floored at Q_floor.
#         AWI > 0: drying down (wetter recent past); AWI < 0: wetting up.
#         The ratio avoids re-encoding today's Q, which is already in the power law.

season_of <- function(d) {
  m <- month(d)
  factor(case_when(m %in% c(12, 1, 2) ~ "Winter (DJF)",
                   m %in% 3:5         ~ "Spring (MAM)",
                   m %in% 6:8         ~ "Summer (JJA)",
                   TRUE               ~ "Fall (SON)"),
         levels = c("Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Fall (SON)"))
}

# Trailing rolling mean over k values, ignoring NA; returns mean and count
roll_mean_na <- function(x, k) {
  v <- cumsum(replace(x, is.na(x), 0))
  n <- cumsum(!is.na(x))
  s <- v - dplyr::lag(v, k, default = 0)
  c <- n - dplyr::lag(n, k, default = 0)
  list(mean = if_else(c > 0, s / c, NA_real_), n = c)
}

drivers <- hydro %>%
  group_by(watershed) %>%
  complete(date = seq(min(date), max(date), by = "day")) %>%   # lags in days
  arrange(date, .by_group = TRUE) %>%
  mutate(Q_prev = lag(Q_mm_day),
         dQ_rel = (Q_mm_day - Q_prev) / Q_prev,
         limb   = case_when(is.na(dQ_rel) | !is.finite(dQ_rel) ~ "steady",
                            dQ_rel >  limb_tol                 ~ "rising",
                            dQ_rel < -limb_tol                 ~ "falling",
                            TRUE                               ~ "steady"),
         Q_ante = roll_mean_na(Q_prev, ante_days)$mean,
         n_ante = roll_mean_na(Q_prev, ante_days)$n,
         Q_ante = if_else(n_ante >= ante_min_days, Q_ante, NA_real_),
         AWI    = log10(pmax(Q_ante, Q_floor) / pmax(Q_mm_day, Q_floor))) %>%
  ungroup() %>%
  select(watershed, date, limb, AWI)

resid_data <- resid_data %>%
  left_join(drivers, by = c("watershed", "date")) %>%
  mutate(limb   = factor(limb, levels = c("rising", "falling", "steady")),
         season = season_of(date))

# 6.1 Table: dispersion and R2 added by each driver -------------------------
# Added R2 when each term joins the power law. All models use the same days
# (those with an AWI value) so values are comparable.
explain_one <- function(d) {
  r2   <- function(f) summary(lm(f, data = d))$r.squared
  base <- r2(lASDN ~ lQ)
  tibble(n_driver_days = nrow(d),
         r2_powerlaw   = base,
         dR2_limb      = r2(lASDN ~ lQ + limb)   - base,
         dR2_season    = r2(lASDN ~ lQ + season) - base,
         dR2_AWI       = r2(lASDN ~ lQ + AWI)    - base,
         dR2_all       = r2(lASDN ~ lQ + limb + season + AWI) - base,
         rho_resid_AWI = cor(d$resid, d$AWI, method = "spearman"))
}

table_dispersion <- dispersion %>%
  left_join(resid_data %>%
              filter(limb != "steady") %>%
              group_by(watershed) %>%
              summarise(resid_rise_minus_fall = mean(resid[limb == "rising"]) -
                                                mean(resid[limb == "falling"]),
                        p_rise_vs_fall = t.test(resid[limb == "rising"],
                                                resid[limb == "falling"])$p.value,
                        .groups = "drop"),
            by = "watershed") %>%
  left_join(resid_data %>%
              filter(!is.na(AWI)) %>%
              group_by(watershed) %>%
              group_modify(~explain_one(.x)) %>%
              ungroup(),
            by = "watershed") %>%
  left_join(resid_data %>%
              group_by(watershed, season) %>%
              summarise(m = mean(resid), .groups = "drop") %>%
              mutate(season = paste0("resid_mean_", str_extract(season, "(?<=\\()[A-Z]+"))) %>%
              pivot_wider(names_from = season, values_from = m),
            by = "watershed") %>%
  mutate(across(where(is.double) & !starts_with("p_"), ~round(.x, 3)),
         p_rise_vs_fall = signif(p_rise_vs_fall, 2))
print(table_dispersion, width = Inf)

write_csv(table_dispersion, file.path(tab_dir, "SI_table_dispersion_drivers.csv"))

# 6.2 Figure: ASDN vs. Q colored by season -------------------------------------
season_cols <- c("Winter (DJF)" = "#0072B2", "Spring (MAM)" = "#009E73",
                 "Summer (JJA)" = "#D55E00", "Fall (SON)"   = "#E69F00")

fig_season <- ggplot(resid_data, aes(x = Q_mm_day, y = ASDN_pct, color = season)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_line(data = pl_lines, aes(x = Q, y = ASDN_hat), color = "black",
            linewidth = 0.8, inherit.aes = FALSE) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_color_manual(values = season_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = "Daily mean discharge (mm/day)", y = "ASDN (% of maximum active length)",
       caption = "Black line: per-site power-law fit (power-law figure). Meteorological seasons.") +
  theme_hydro +
  theme(legend.position = "bottom")
fig_season

ggsave(file.path(fig_dir, "SI_ASDN_Q_by_season.png"), fig_season,
       width = 10, height = 9, dpi = 300)

# 6.3 Figure: residuals vs. Q colored by season --------------------------------
# Same residuals as the residuals-vs-Q figure, colored by season. Seasonal
# mean residuals are in the dispersion table (resid_mean_* columns).
fig_resid_season <- ggplot(resid_data, aes(x = Q_mm_day, y = resid, color = season)) +
  geom_rect(data = dispersion, aes(xmin = 0, xmax = Inf, ymin = -resid_sd, ymax = resid_sd),
            fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 1, alpha = 0.5) +
  geom_text(data = table_dispersion,
            aes(x = 0, y = Inf, label = sprintf("ΔR² season = %.2f", dR2_season)),
            hjust = -0.1, vjust = 1.4, size = 3, color = "black", inherit.aes = FALSE) +
  scale_x_log10(labels = log_labels) +
  scale_color_manual(values = season_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = "Daily mean discharge (mm/day)",
       y = "Residual, log10(ASDN observed / ASDN predicted)",
       caption = paste("Residuals from the per-site power-law fits. Gray band: ±1 SD.",
                       "ΔR² season: R² added when season joins the power law.",
                       "Meteorological seasons.")) +
  theme_hydro +
  theme(legend.position = "bottom")
fig_resid_season

ggsave(file.path(fig_dir, "SI_residuals_vs_Q_by_season.png"), fig_resid_season,
       width = 10, height = 9, dpi = 300)

# 6.4 Figure: residuals vs. Q colored by day of year ---------------------------
# Cyclic color scale so Dec 31 and Jan 1 match. Season colors are anchored near
# each season's center (winter at the year boundary).
fig_resid_doy <- resid_data %>%
  mutate(doy = yday(date)) %>%
  ggplot(aes(x = Q_mm_day, y = resid, color = doy)) +
  geom_rect(data = dispersion, aes(xmin = 0, xmax = Inf, ymin = -resid_sd, ymax = resid_sd),
            fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 1, alpha = 0.6) +
  scale_x_log10(labels = log_labels) +
  scale_color_gradientn(colors = unname(season_cols[c(1, 2, 3, 4, 1)]),
                        values = scales::rescale(c(1, 105, 196, 288, 366)),
                        limits = c(1, 366),
                        breaks = yday(as_date(paste0("2021-", c("01", "04", "07", "10"), "-01"))),
                        labels = c("Jan", "Apr", "Jul", "Oct"),
                        name = "Day of year") +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = "Daily mean discharge (mm/day)",
       y = "Residual, log10(ASDN observed / ASDN predicted)",
       caption = "Residuals from the per-site power-law fits. Gray band: ±1 SD. Cyclic day-of-year color scale.") +
  theme_hydro +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"))
fig_resid_doy

ggsave(file.path(fig_dir, "SI_residuals_vs_Q_by_doy.png"), fig_resid_doy,
       width = 10, height = 9, dpi = 300)

# 6.5 Figure: residuals vs. antecedent wetness ---------------------------------
fig_awi <- resid_data %>%
  filter(!is.na(AWI)) %>%
  ggplot(aes(x = AWI, y = resid, color = watershed)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 0.9, alpha = 0.35) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.7) +
  geom_text(data = table_dispersion,
            aes(x = -Inf, y = Inf,
                label = sprintf("ρ = %.2f\nΔR² = %.2f", rho_resid_AWI, dR2_AWI)),
            hjust = -0.1, vjust = 1.2, size = 3, color = "black", inherit.aes = FALSE) +
  scale_color_manual(values = col_line, guide = "none") +
  facet_wrap(~watershed, ncol = 3, labeller = site_lab) +
  labs(x = paste0("Antecedent wetness index, log10(mean Q previous ", ante_days, " d / Q today)"),
       y = "Residual, log10(ASDN observed / predicted)",
       caption = paste("AWI > 0: drying down; AWI < 0: wetting up.",
                       "ρ: Spearman correlation of residual with AWI.",
                       "ΔR²: R² added when AWI joins the power law.")) +
  theme_hydro
fig_awi

ggsave(file.path(fig_dir, "SI_residuals_vs_AWI.png"), fig_awi,
       width = 10, height = 9, dpi = 300)

# 7.0 Main figure 1: elasticity and dispersion at an example watershed ---------
# (a) Q and ASDN time series, (b) log-log ASDN vs. Q with the power law
# (elasticity), (c) residuals vs. Q (dispersion), (d, e) STIC wet/dry maps on two
# days with similar Q but opposite residuals: same outlet flow, different network.
#
# Node states come from the raw STIC files (R/spatial_functions.R), because the
# merged file only has network totals. Raw daily wet counts track streamDAG
# wet_network_nodes closely (r = 0.97 at WHR) but match exactly on only ~1/3 of
# days (streamDAG has all 20 nodes on every day and starts before the raw STIC
# record, so it likely gap-fills). Map days are therefore days where the two
# counts match exactly (checked below).
library(sf)
source("R/spatial_functions.R")

ex_ws    <- "WHR"
map_days <- as_date(c("2024-05-16", "2024-06-15"))
# Days chosen with all 20 STICs reporting, wet count = streamDAG, ASDN < 100%,
# Q within ~0.15 log units, residuals of opposite sign. Spring vs. early summer
# 2024: a typical seasonal (ET-driven) contrast. Stronger alternative:
# 2023-05-15 (Q 0.077 mm/d, ASDN 92%) vs. 2024-07-17 (Q 0.090 mm/d, ASDN 17%),
# but the latter is the first day of a storm after 6 days of zero flow.

# 7.1 Spatial data -------------------------------------------------------------
ex_shed   <- read_zip_shape(ex_ws, "watershed")
ex_crs    <- st_crs(ex_shed)
ex_net    <- read_zip_shape(ex_ws, "streamnetwork") %>% st_set_crs(ex_crs)  # same UTM 16N, no EPSG tag
ex_stic   <- read_stic_sites(ex_ws, ex_crs)
ex_outlet <- read_envi_sites(ex_ws, ex_crs) %>% filter(distance_from_outlet == 0)
ex_daily  <- read_stic_daily(ex_ws)

# Q, ASDN, residual, and node counts on the map days
map_info <- fit_data %>%
  filter(watershed == ex_ws, date %in% map_days) %>%
  select(date, Q_mm_day, ASDN_pct, resid) %>%
  left_join(read_csv(merge_file, show_col_types = FALSE) %>%
              filter(watershed == ex_ws) %>%
              select(date, wet_nodes_streamDAG = wet_network_nodes),
            by = "date") %>%
  left_join(ex_daily %>%
              group_by(date) %>%
              summarise(n_stic = n(), wet_nodes_stic = sum(wet), .groups = "drop"),
            by = "date") %>%
  arrange(date) %>%
  mutate(day_lab = c("i", "ii")[seq_len(n())],
         title   = sprintf("(%s) %s: Q = %.2f mm/day, ASDN = %.0f%%",
                           c("d", "e")[seq_len(n())], format(date, "%d %b %Y"),
                           Q_mm_day, ASDN_pct))
print(map_info, width = Inf)
if (any(map_info$wet_nodes_stic != map_info$wet_nodes_streamDAG)) {
  warning("Raw STIC wet count differs from streamDAG on a map day")
}

# 7.2 Panels -------------------------------------------------------------------
ex_col  <- col_line[[ex_ws]]
ex_pts  <- filter(fit_data, watershed == ex_ws, use_in_fit)
ex_fit  <- filter(pl_fits, watershed == ex_ws)
ex_disp <- filter(dispersion, watershed == ex_ws)
ex_dates <- range(filter(hydro, watershed == ex_ws)$date)

main_theme <- theme_hydro +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "plain"))

# Highlighted map days on the scatter panels
mark_days <- function(y) {
  list(geom_point(data = map_info, aes(x = Q_mm_day, y = {{ y }}),
                  shape = 21, fill = "white", color = "black", size = 3, stroke = 0.9),
       geom_text(data = map_info, aes(x = Q_mm_day, y = {{ y }}, label = day_lab),
                 nudge_x = 0.13, hjust = 0, vjust = 0.5, size = 3.3, fontface = "bold"))
}

# (a) Time series
p1a <- plot_ts_site(ex_ws, "", show_x = TRUE, show_ylab = TRUE) +
  geom_vline(xintercept = map_info$date, linetype = "dashed", linewidth = 0.4) +
  annotate("label", x = map_info$date, y = 100, label = map_info$day_lab,
           size = 3, fontface = "bold", label.size = 0, label.padding = unit(0.1, "lines")) +
  annotate("text", x = ex_dates[1], y = 100, label = "(a)", hjust = 0.2, vjust = 1, size = 3.5) +
  scale_x_date(limits = ex_dates, date_breaks = "6 months", date_labels = "%b %Y",
               expand = expansion(mult = 0.01)) +
  labs(title = NULL, x = NULL) +
  main_theme +
  theme(axis.title.y.left  = element_text(color = "black"),
        axis.title.y.right = element_text(color = "grey25"))

# (b) Elasticity
p1b <- ggplot(ex_pts, aes(x = Q_mm_day, y = ASDN_pct)) +
  geom_point(shape = 21, fill = ex_col, color = "grey40", stroke = 0.2, size = 1.4, alpha = 0.8) +
  geom_line(data = filter(pl_lines, watershed == ex_ws), aes(x = Q, y = ASDN_hat), linewidth = 0.9) +
  mark_days(ASDN_pct) +
  annotate("text", x = Inf, y = 0, hjust = 1.05, vjust = -0.5, size = 3,
           label = sprintf("Elasticity b = %.2f\nR² = %.2f", ex_fit$elasticity, ex_fit$r2)) +
  annotate("text", x = 0, y = Inf, label = "(b)", hjust = -0.3, vjust = 1.4, size = 3.5) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10(breaks = c(5, 10, 20, 50, 100)) +
  labs(x = "Daily mean discharge (mm/day)", y = "ASDN (%)") +
  main_theme

# (c) Dispersion
p1c <- ggplot(ex_pts, aes(x = Q_mm_day, y = resid)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -ex_disp$resid_sd, ymax = ex_disp$resid_sd,
           fill = "grey80", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(shape = 21, fill = ex_col, color = "grey40", stroke = 0.2, size = 1.4, alpha = 0.8) +
  mark_days(resid) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, size = 3,
           label = sprintf("Dispersion (SD) = %.2f\n(×%.2f)", ex_disp$resid_sd, ex_disp$factor_sd)) +
  annotate("text", x = 0, y = Inf, label = "(c)", hjust = -0.3, vjust = 1.4, size = 3.5) +
  scale_x_log10(labels = log_labels) +
  labs(x = "Daily mean discharge (mm/day)", y = "Residual, log10(obs/pred)") +
  main_theme

# (d, e) Network maps
# Reach states between STICs (rules in classify_reaches(), R/spatial_functions.R):
# a reach is wet if the STICs on both ends are wet; at a confluence the shared
# reach is wet if the downstream STIC and either upstream STIC are wet; headwater
# reaches follow the STIC below them and the outlet reach the STIC above it.
ex_graph <- build_stream_graph(ex_net, ex_stic, ex_outlet)

plot_network_day <- function(day, title, scale_bar = FALSE) {
  wet_by_site <- ex_daily %>% filter(date == day) %>% select(siteId, wet) %>% deframe()
  reaches <- classify_reaches(ex_graph, wet_by_site) %>%
    filter(!is.na(wet)) %>%   # drops short unmonitored stubs at the outlet (no STIC above or below)
    mutate(state = factor(if_else(wet, "Wet", "Dry"), levels = c("Wet", "Dry")))
  bb <- st_bbox(ex_shed)
  p <- ggplot() +
    geom_sf(data = ex_shed, fill = "grey96", color = "grey40", linewidth = 0.4) +
    geom_sf(data = reaches, aes(color = state), linewidth = 1.2, lineend = "round") +
    geom_sf(data = ex_outlet, shape = 24, fill = "black", size = 2.2) +
    scale_color_manual(values = c(Wet = "#2166AC", Dry = "#F4A582"), drop = FALSE,
                       name = "Reach state", na.value = "grey70") +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(size = 9.5, margin = margin(b = 3)),
          plot.margin = margin(4, 4, 4, 4))
  if (scale_bar) {
    x0 <- bb[["xmin"]] + 30; y0 <- bb[["ymin"]] + 30
    p <- p +
      annotate("segment", x = x0, xend = x0 + 250, y = y0, yend = y0, linewidth = 0.8) +
      annotate("text", x = x0 + 125, y = y0 + 25, label = "250 m", size = 2.6)
  }
  p
}

p1d <- plot_network_day(map_info$date[1], map_info$title[1], scale_bar = TRUE)
p1e <- plot_network_day(map_info$date[2], map_info$title[2])

# 7.3 Assemble -----------------------------------------------------------------
fig_main_1 <- p1a / (p1b | p1c) / (p1d | p1e) +
  plot_layout(heights = c(1, 1.15, 0.95), guides = "collect") &
  theme(legend.position = "bottom")
fig_main_1

ggsave(file.path(fig_dir, "main_WHR_elasticity_dispersion.png"), fig_main_1,
       width = 7.5, height = 9, dpi = 300)

# 8.0 Main figure 2: elasticity and dispersion across sites --------------------
# (a) Daily ASDN vs. Q for all sites with each site's fitted power law
# (b) Elasticity vs. dispersion, one point per site (bars: 95% CI on b)
# Note: the CIs on b assume independent days; with autocorrelated daily data
# they are too narrow and are shown only as a relative guide.
library(ggrepel)

region_shapes <- c("Great Plains" = 21, "Mountain West" = 22, "Southeast" = 24)

site_metrics <- pl_fits %>%
  select(watershed, elasticity, b_ci_low, b_ci_high, r2) %>%
  left_join(dispersion, by = "watershed") %>%
  left_join(site_info %>% mutate(watershed = factor(watershed, levels = site_order)),
            by = "watershed")

# (a) Daily points (faint) under the fitted lines (white halo for contrast)
pl_ends <- pl_lines %>% group_by(watershed) %>% slice_max(lQ, n = 1) %>% ungroup()

p2a <- fit_data %>%
  filter(use_in_fit) %>%
  ggplot(aes(x = Q_mm_day, y = ASDN_pct, color = watershed)) +
  geom_point(size = 0.7, alpha = 0.35, shape = 16) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(data = pl_lines, aes(x = Q, y = ASDN_hat, group = watershed),
            color = "white", linewidth = 2.4) +
  geom_line(data = pl_lines, aes(x = Q, y = ASDN_hat), linewidth = 1.1) +
  geom_text_repel(data = pl_ends, aes(x = Q, y = ASDN_hat, label = watershed),
                  size = 3, fontface = "bold", direction = "y", hjust = 0, nudge_x = 0.25,
                  segment.color = "grey60", segment.size = 0.3, min.segment.length = 0,
                  bg.color = "white", bg.r = 0.12, seed = 1) +
  annotate("text", x = 0, y = Inf, label = "(a)", hjust = -0.3, vjust = 1.4, size = 3.5) +
  scale_x_log10(labels = log_labels, expand = expansion(mult = c(0.02, 0.17))) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_color_manual(values = col_line, guide = "none") +
  labs(x = "Daily mean discharge (mm/day)", y = "ASDN (% of maximum active length)") +
  main_theme

# (b) Elasticity vs. dispersion
p2b <- ggplot(site_metrics, aes(x = resid_sd, y = elasticity)) +
  geom_errorbar(aes(ymin = b_ci_low, ymax = b_ci_high), width = 0, color = "grey40") +
  geom_point(aes(fill = watershed, shape = region), size = 3.4, color = "black", stroke = 0.4) +
  geom_text_repel(aes(label = watershed), size = 3, fontface = "bold", seed = 1,
                  point.padding = unit(0.4, "lines"), min.segment.length = 1) +
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.3, vjust = 1.4, size = 3.5) +
  scale_fill_manual(values = col_line, guide = "none") +
  scale_shape_manual(values = region_shapes, name = NULL) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  guides(shape = guide_legend(override.aes = list(fill = "grey70", size = 3))) +
  labs(x = "Dispersion (SD of log10 residuals)", y = "Elasticity (b)") +
  main_theme +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03),
        legend.justification = c(1, 0), legend.background = element_rect(fill = NA),
        legend.key.size = unit(0.9, "lines"), legend.text = element_text(size = 8))

# 8.1 Assemble -----------------------------------------------------------------
fig_main_2 <- p2a | p2b
fig_main_2

ggsave(file.path(fig_dir, "main_cross_site_elasticity_dispersion.png"), fig_main_2,
       width = 7.5, height = 3.9, dpi = 300)
