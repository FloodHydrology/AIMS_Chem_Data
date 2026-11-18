#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Title: Spatial helper functions for AIMS watersheds
# Date: 11/18/2026
# Coder: Nate Jones (cnjones7@ua.edu)
# Purpose: Read watershed boundaries, stream networks, STIC locations, and daily
#          STIC wet/dry states for one watershed, keyed by its code (e.g., "WHR").
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Expected layout (HydroShare download, R/I_Download.R):
#   <data_dir>/<ws>_ENVI/watershed_<ws>.zip       watershed boundary shapefile
#   <data_dir>/<ws>_ENVI/streamnetwork_<ws>.zip   stream network shapefile
#   <data_dir>/<ws>_ENVI/ENVI_*_<ws>.xlsx         site attributes ("Final Data")
#   <data_dir>/<ws>_STIC/STIC_*.csv               15-min STIC records (wetdry 0/1)
#   <data_dir>/<ws>_STIC/*ReadME.xlsx             STIC coordinates ("Site_Info")
# Requires: tidyverse, lubridate, readxl, sf
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read a zipped shapefile ("watershed" or "streamnetwork") for one watershed
read_zip_shape <- function(ws, layer, data_dir = "aims_data") {
  zip_file <- file.path(data_dir, paste0(ws, "_ENVI"), paste0(layer, "_", ws, ".zip"))
  if (!file.exists(zip_file)) stop("Missing ", zip_file)
  out_dir <- file.path(tempdir(), paste0(layer, "_", ws))
  unzip(zip_file, exdir = out_dir)
  shp <- list.files(out_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  st_read(shp[1], quiet = TRUE)
}

# STIC locations (Site_Info sheet of the STIC ReadMe), projected to crs
read_stic_sites <- function(ws, crs, data_dir = "aims_data") {
  readme <- list.files(file.path(data_dir, paste0(ws, "_STIC")),
                       pattern = "ReadME\\.xlsx$", full.names = TRUE)
  if (length(readme) == 0) stop("No STIC ReadMe found for ", ws)
  read_excel(readme[1], sheet = "Site_Info") %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(crs)
}

# All environmental sites (ENVI "Final Data" sheet), projected to crs.
# The outlet is the site with distance_from_outlet == 0.
read_envi_sites <- function(ws, crs, data_dir = "aims_data") {
  envi <- list.files(file.path(data_dir, paste0(ws, "_ENVI")),
                     pattern = "^ENVI_.*\\.xlsx$", full.names = TRUE)
  if (length(envi) == 0) stop("No ENVI file found for ", ws)
  read_excel(envi[1], sheet = "Final Data") %>%
    mutate(distance_from_outlet = as.numeric(distance_from_outlet)) %>%
    st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
    st_transform(crs)
}

# Daily STIC wet/dry state per site from the 15-min records.
#   frac_wet = fraction of readings that day classified wet (wetdry == 1)
#   wet      = frac_wet >= wet_threshold
# Days are UTC calendar days (datetime is stored in UTC).
read_stic_daily <- function(ws, data_dir = "aims_data", wet_threshold = 0.5) {
  files <- list.files(file.path(data_dir, paste0(ws, "_STIC")),
                      pattern = "^STIC_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No STIC data files found for ", ws)
  map_dfr(files, ~read_csv(.x, col_types = cols_only(datetime = "c", siteId = "c",
                                                      wetdry = "i"))) %>%
    filter(!is.na(wetdry)) %>%
    mutate(date = as_date(ymd_hms(datetime))) %>%
    group_by(siteId, date) %>%
    summarise(n_obs = n(), frac_wet = mean(wetdry), .groups = "drop") %>%
    mutate(wet = frac_wet >= wet_threshold)
}

# Directed stream graph with STICs attached, for classifying reaches wet/dry.
# The network is split at every vertex; vertices shared by lines (junctions)
# join the graph. Flow direction is toward the vertex nearest the outlet.
# For each small edge:
#   D = nearest STIC downstream (first STIC met walking toward the outlet)
#   U = nearest STIC(s) upstream (first STIC met on each upstream path)
# Returns an sf of edges (LINESTRING) with list-column U and column D.
# Requires igraph (called with igraph::, not attached).
build_stream_graph <- function(net, sites, outlet, site_col = "siteId") {
  co <- st_coordinates(st_cast(st_geometry(net), "LINESTRING")) %>%
    as_tibble() %>%
    mutate(key = paste(round(X, 2), round(Y, 2)))
  v <- co %>% distinct(key, .keep_all = TRUE) %>% transmute(key, X, Y, vid = row_number())
  e <- co %>%
    left_join(select(v, key, vid), by = "key") %>%
    group_by(L1) %>%
    transmute(a = vid, b = lead(vid)) %>%
    ungroup() %>%
    filter(!is.na(b), a != b) %>%
    distinct(a, b) %>%
    mutate(len = sqrt((v$X[a] - v$X[b])^2 + (v$Y[a] - v$Y[b])^2))
  g <- igraph::graph_from_data_frame(transmute(e, from = a, to = b, weight = len),
                                     directed = FALSE, vertices = select(v, vid))
  if (igraph::components(g)$no != 1) warning("Stream network is not fully connected")

  # Snap outlet and sites to their nearest vertex
  v_sf <- st_as_sf(v, coords = c("X", "Y"), crs = st_crs(net))
  out_vid <- st_nearest_feature(outlet, v_sf)[1]
  v$site  <- NA_character_
  v$site[st_nearest_feature(sites, v_sf)] <- sites[[site_col]]

  # Distance to outlet orients each edge: child (upstream) -> parent (downstream)
  v$dist <- as.numeric(igraph::distances(g, v = out_vid)[1, ])
  e <- e %>%
    mutate(child  = if_else(v$dist[a] > v$dist[b], a, b),
           parent = if_else(v$dist[a] > v$dist[b], b, a))
  parent_of <- rep(NA_integer_, nrow(v)); parent_of[e$child] <- e$parent
  children  <- split(e$child, factor(e$parent, levels = seq_len(nrow(v))))

  # Nearest STIC downstream of each vertex (outlet-first order)
  D <- rep(NA_character_, nrow(v))
  for (i in order(v$dist)) {
    D[i] <- if (!is.na(v$site[i])) v$site[i] else if (is.na(parent_of[i])) NA_character_ else D[parent_of[i]]
  }
  # Nearest STIC(s) upstream of each vertex (headwater-first order)
  U <- vector("list", nrow(v))
  for (i in order(v$dist, decreasing = TRUE)) {
    U[[i]] <- if (!is.na(v$site[i])) v$site[i] else unique(unlist(U[children[[i]]]))
  }

  geom <- st_sfc(lapply(seq_len(nrow(e)), function(k) {
    st_linestring(rbind(c(v$X[e$child[k]],  v$Y[e$child[k]]),
                        c(v$X[e$parent[k]], v$Y[e$parent[k]])))
  }), crs = st_crs(net))
  st_sf(len = e$len, U = I(U[e$child]), D = D[e$parent], geometry = geom)
}

# Wet/dry state of each edge from one day's STIC states (named logical vector)
#   both neighbors known: wet if D is wet and any upstream STIC is wet
#   headwater (no U):     wet if D is wet
#   outlet reach (no D):  wet if any upstream STIC is wet
classify_reaches <- function(graph_edges, wet_by_site) {
  graph_edges$wet <- map2_lgl(graph_edges$U, graph_edges$D, function(u, d) {
    up <- if (length(u)) any(wet_by_site[u]) else NA
    dn <- if (!is.na(d)) wet_by_site[[d]] else NA
    if (length(u) == 0) dn else if (is.na(d)) up else dn & up
  })
  graph_edges
}
