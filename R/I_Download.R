# =============================================================================
# Download AIMS data from HydroShare  (v2 — per-file API)
#
# For each resource: query /hsapi/resource/{id}/files/ to get the file list,
# then download each file directly into aims_data/<LABEL>/.
# Files are .xlsx / .zip (shapefiles) / .csv depending on resource.
# All resources are public (CC-BY, Published). No auth needed.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
for (p in c("httr", "jsonlite")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(httr); library(jsonlite)

# ---- Output directory -------------------------------------------------------
out_dir <- "aims_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Resource manifest (label -> HydroShare resource id) --------------------
resources <- read.csv(text = "
label,id
TAL_ENVI,81c003a7b8474d63a31641a4f375fd18
PRF_ENVI,656211b1a1484433a3bc524fb968b4bd
WHR_ENVI,126d2c7b1c8d4889a8ccc454d387b0d8
TAL_METS,281cd7627629481dbdc7d4ccf6fcfbcc
PRF_METS,4089918c0a494bfeb19be0421a33d297
WHR_METS,33823d8603ce439fbba48fbcbba22da4
TAL_STIC,ff306bec9fb24e52aa809dbb4d074731
PRF_STIC,d57338ebfb0240f58e8de37ddacf9426
WHR_STIC,dc623510ed1847f8abe1275904472c44
TAL_PRES,93e2861410e647d9a710eea036832dbe
PRF_PRES,a45b5e24dafc4a76a665405664afada7
WHR_PRES,bc34c8b51c514bf4a6e0a44493bf8ca3
TAL_DISC,fc7ae2d28e3c481d805902a79af90a95
PRF_DISC,043fc07f0c3b47bcabbd0bf5600d929f
WHR_DISC,535797126b134ceaab9838df0ca00885
TAL_EXOS,9e47c3dfc549446e80173dfe6ac48365
PRF_EXOS,aeb715ff9c7b4f1098bdebc0fd9e9551
WHR_EXOS,a7bdb79e06684db2886a257ec614018a
TAL_SCAN,cea7ec0e055f49ef9f55fc61caffc52a
TAL_YSIS,e36dc69dca0e4fbc969e7ae6137f3744
PRF_YSIS,7fb2e1872cb840bcbcd8bd4e1ea12185
WHR_YSIS,a9394fd2e0d748fbb3ca5c36b451c15f
TAL_DISL,0e7ad0451bdc45d2b0a51bb538a10909
PRF_DISL,d52b989e537349019842dba236627b66
WHR_DISL,eedcfcb232ee45a6915bd26c68e301e8
TAL_WAIS,5fffa420810f4945a7ee8d3f8bda3ad2
PRF_WAIS,2088b266609c4ac58195a9390598633e
WHR_WAIS,3ac22e1bc2c547e8a8af4ebd753323e4
TAL_TSSS,b230a7995b06498cacd28106a3be0f35
PRF_TSSS,3eaacf0102594482ae2451c60745d7e6
WHR_TSSS,1284a362f1f9410b87d91598afc53c83
TAL_NUTR,730b486c0ef14d78b678963ffecc1a39
PRF_NUTR,165f2b4d1903485d82304bbc55ecd715
WHR_NUTR,c0008581efb741dda4156fa887c16eb5
TAL_DOCS,e80e4db42de940aa9fe18667dddebec4
PRF_DOCS,1efa655d91fe43c58f8acbf0f52545c8
WHR_DOCS,8b750838affc438e88bcb2cd0dfd5dbf
TAL_DOMS,aa792fa579f5443bba4376008da9f48e
PRF_DOMS,da3766f455944ef0a8613c20d2870d38
WHR_DOMS,b6142a7bfabe4be988733e2c59bd8533
TAL_ANIO,0decb1efb3a34e88b39b64dbb6369743
PRF_ANIO,e025fb27b18141beab4cebda71528efc
WHR_ANIO,3a783ac086a74e6987dfefb870fb8cb3
TAL_CAIO,dc0434b19c834941aa56449af0f6ce9b
PRF_CAIO,0495c3eabb474b1190211fd278b4b467
WHR_CAIO,eb3f2e78492f4ec9bd5a11791712a6f9
TAL_MIMS,5ff9056710d04917bd6891b46496d7b0
PRF_MIMS,036b5916526347bc8bad0ad61559fb9e
WHR_MIMS,c8b4f7ebda48424fad3d709a1b9372aa
TAL_GHGS,34b55fc99e94410f8db6766511b448bb
SE_MIME,3161225427d8472d9f347068e1afab61
TAL_MAME,549b107d949e43cba49adadfdc9b0c15
PRF_MAME,21421686430f42ca9e9936fac26fffd3
WHR_MAME,b2f68f520074419f8e556585daa5b371
SE_AFDM,df5dff9fd883414a8bf91ddeb268e514
SE_CHLA,cd2852e4a0ca4e8d8d65dd3bcd7bd8ad
SE_MACR,8f8d336d073343e7af1197d1ce6b6085
TAL_EEAS_A2,433c5de6768d4ad89f0027ad2101dcda
TAL_EEAS_A3,eb7624d386584c1fb5468ed376487552
PRF_EEAS,3b2886a7bade49dabc5a7d1413b73681
WHR_EEAS,b4dafa88679a444da26261d4c47ee784
", stringsAsFactors = FALSE, strip.white = TRUE)
resources <- resources[nzchar(resources$label), ]

cat(sprintf("Manifest: %d resources.\n\n", nrow(resources)))

# ---- Download one resource --------------------------------------------------
download_one <- function(label, id, out_dir, timeout_sec = 600) {
  target <- file.path(out_dir, label)
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Ask the API what files this resource holds.
  list_url <- sprintf("https://www.hydroshare.org/hsapi/resource/%s/files/", id)
  lr <- tryCatch(GET(list_url, timeout(120)), error = function(e) e)
  if (inherits(lr, "error") || status_code(lr) != 200) {
    cat(sprintf("  [%s] LIST FAILED (%s)\n", label,
                if (inherits(lr, "error")) conditionMessage(lr)
                else paste("HTTP", status_code(lr))))
    return(data.frame(label = label, n_files = 0, ok = FALSE))
  }
  
  parsed <- fromJSON(rawToChar(lr$content), simplifyDataFrame = TRUE)
  files  <- parsed$results
  if (is.null(files) || nrow(files) == 0) {
    cat(sprintf("  [%s] no files listed\n", label))
    return(data.frame(label = label, n_files = 0, ok = FALSE))
  }
  
  # 2. Download each file. Force https to avoid redirect hiccups.
  n_ok <- 0
  for (j in seq_len(nrow(files))) {
    fname <- files$file_name[j]
    furl  <- sub("^http://", "https://", files$url[j])
    dst   <- file.path(target, fname)
    dr <- tryCatch(
      GET(furl, write_disk(dst, overwrite = TRUE), timeout(timeout_sec)),
      error = function(e) e
    )
    if (!inherits(dr, "error") && status_code(dr) == 200) {
      n_ok <- n_ok + 1
    } else {
      cat(sprintf("    ! %s/%s failed\n", label, fname))
    }
  }
  
  cat(sprintf("  [%s] %d/%d files\n", label, n_ok, nrow(files)))
  data.frame(label = label, n_files = nrow(files), ok = (n_ok == nrow(files)))
}

# ---- Run --------------------------------------------------------------------
res_list <- vector("list", nrow(resources))
for (i in seq_len(nrow(resources))) {
  res_list[[i]] <- download_one(resources$label[i], resources$id[i], out_dir)
  Sys.sleep(0.5)  # be polite
}
summary_df <- do.call(rbind, res_list)

# ---- Summary ----------------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("Fully succeeded: %d / %d resources\n",
            sum(summary_df$ok), nrow(summary_df)))
cat(sprintf("Total files downloaded: %d\n", sum(summary_df$n_files[summary_df$ok])))
bad <- summary_df$label[!summary_df$ok]
if (length(bad)) {
  cat("\nProblem resources (re-run these):\n")
  cat(paste0("  ", bad, collapse = "\n"), "\n")
}
cat(sprintf("\nData is in: %s/<LABEL>/\n", normalizePath(out_dir)))