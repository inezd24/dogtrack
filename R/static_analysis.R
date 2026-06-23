# =============================================================================
# R/static_analysis.R
# Package : dogtrack
# Purpose : Load static GPS test files and derive site-specific thresholds.
#
#   Static tests = GPS collar placed motionless on the ground for an extended
#   period. Any positional variation is pure device + environmental noise.
#
#   Derived thresholds (per field site):
#     dist_noise_m     -- p95 Haversine distance from centroid (metres)
#                        used as the distance gate in filter_angle()
#     height_threshold -- p99.9 ellipsoidal HEIGHT + height_buffer_m
#                        used as the ceiling in filter_height()
#
# Exported functions:
#   load_static_tests()      -- reads static test files, tags by site and device
#   derive_site_thresholds() -- computes dist_noise_m + height_threshold per site
#
# NOTE: load_static_tests() depends on internal helpers defined in
# R/load_data.R: .detect_device(), .read_one_columbus(), .read_one_catlogger(),
# .read_one_other(). Both files must be part of the same package build.
# =============================================================================


#' Load static GPS test files
#'
#' Scans a directory for static test recordings -- stationary GPS logs used
#' to derive site-specific accuracy and precision thresholds -- and combines
#' them into a single tidy tibble. Supports Columbus, CatLog, and other
#' devices via the same `device_type` mechanism as [load_dog_data()].
#'
#' @param data_dir Character. Path to the directory containing the GPS data
#'   files (dog CSVs and static test CSVs may be mixed or in sub-folders).
#' @param prefix_map A tibble as returned by [build_prefix_map()], with
#'   columns `prefix_2` and `field_site`.
#' @param file_pattern Character (regex). Identifies which files in
#'   `data_dir` are static test files. Matched case-insensitively against the
#'   file's basename. Defaults to `"^static"`, matching the convention used
#'   so far (e.g. `"staticUM155.CSV"`, `"Static_TB03_columbus.csv"`). Override
#'   this if your static files use a different naming convention or live
#'   under a differently-named subfolder pattern.
#' @param device_type Character. One of `"columbus"`, `"catlogger"`,
#'   `"other"`, or `"mix"` -- see [load_dog_data()] for full details on each.
#' @param col_aliases Named list, or `NULL`. Required when `device_type =
#'   "other"`; see [load_dog_data()].
#'
#' @return A tibble with one row per GPS fix from static test recordings,
#'   containing:
#'   - `static_test_id` : character -- filename (without extension) up to
#'                        and including the site prefix and number (e.g.
#'                        `"Static_TB03"`); any trailing device label is
#'                        excluded from this id and captured separately in
#'                        `device`. NOTE: this is no longer guaranteed unique
#'                        per file -- if the same prefix+number combination
#'                        appears with different device labels (e.g.
#'                        `"Static_TB03_columbus"` and
#'                        `"Static_TB03_catlogger"`), both share this
#'                        `static_test_id`. Use `session_id` (below) when you
#'                        need a per-file, per-device grouping key.
#'   - `session_id`     : character -- the full original filename (without
#'                        extension), guaranteed unique per file. Used by
#'                        [derive_site_thresholds()] for per-session centroid
#'                        computation, so that two device recordings of the
#'                        same physical static test location are kept as
#'                        separate sessions rather than pooled together.
#'   - `device`         : character, or `NA` -- trailing text from the
#'                        filename after the prefix and number (e.g.
#'                        `"columbus"` from `"Static_TB03_columbus"`); `NA`
#'                        if no trailing text is present
#'   - `device_type`    : character -- `"columbus"`, `"catlogger"`, or
#'                        `"other"`, indicating which reader was used to
#'                        parse the file (distinct from `device`, which is
#'                        purely filename-derived text and may or may not
#'                        agree with this)
#'   - `field_site`     : character -- looked up from `prefix_map` via the
#'                        2-letter prefix found in the filename
#'   - `prefix_2`       : character -- the 2-letter site prefix extracted
#'                        from the filename
#'   - `datetime`       : POSIXct -- parsed from DATE + TIME columns
#'   - `lat`            : numeric -- decimal degrees
#'   - `lon`            : numeric -- decimal degrees
#'   - `HEIGHT`         : numeric -- ellipsoidal height (metres, WGS84)
#'   - `HDOP`           : numeric
#'
#'   If no static files are found, a zero-row tibble with these columns is
#'   returned (with a warning), so downstream code can proceed without
#'   special-casing an empty result.
#'
#' @details
#' **File discovery** matches both `.CSV` and `.csv` extensions. Which files
#' count as "static" is controlled by `file_pattern` (default `"^static"`,
#' case-insensitive), checked against the file's basename.
#'
#' **Filename parsing** expects a 2-letter site prefix matching
#' `prefix_map$prefix_2` immediately followed by a number (one or more
#' digits, e.g. `"7"` or `"03"` -- any digit count is accepted). This
#' prefix+number combination becomes `static_test_id`. Anything trailing
#' after that (separated by `_`, `-`, or directly adjacent) is treated as a
#' device label and stored in `device`; if nothing trails, `device` is `NA`.
#' Examples:
#' \itemize{
#'   \item `"staticUA7.CSV"` -> `static_test_id = "staticUA7"`,
#'     `prefix_2 = "UA"`, `device = NA`
#'   \item `"Static_TN30.csv"` -> `static_test_id = "Static_TN30"`,
#'     `prefix_2 = "TN"`, `device = NA`
#'   \item `"Static_TB03_columbus.CSV"` -> `static_test_id = "Static_TB03"`,
#'     `prefix_2 = "TB"`, `device = "columbus"`
#' }
#' Files where no 2-letter prefix can be matched against `prefix_map`, or
#' where the prefix matches but isn't found in `prefix_map`, are skipped
#' with a warning.
#'
#' **Device reading** uses the same `device_type` mechanism as
#' [load_dog_data()]: `"columbus"`/`"catlogger"`/`"other"` apply one reader
#' to every file, `"mix"` auto-detects per file (Columbus vs CatLog only --
#' `"other"`-type static files must be loaded in a separate call). See
#' [load_dog_data()] for the per-device column requirements and behaviour.
#'
#' @examples
#' \dontrun{
#' meta     <- load_metadata("data/metadata.xlsx")
#' pfx_map  <- build_prefix_map(meta)
#'
#' # Default convention, all files Columbus
#' statics <- load_static_tests("data/static/", pfx_map, device_type = "columbus")
#'
#' # Files live in a differently-named folder/convention
#' statics <- load_static_tests(
#'   "data/calibration/", pfx_map,
#'   file_pattern = "^calib", device_type = "mix"
#' )
#' }
#'
#' @export
load_static_tests <- function(data_dir, prefix_map,
                              file_pattern = "^static",
                              device_type,
                              col_aliases = NULL) {
  
  valid_device_types <- c("columbus", "catlogger", "other", "mix")
  if (!device_type %in% valid_device_types) {
    stop(
      "load_static_tests(): `device_type` must be one of: ",
      paste(valid_device_types, collapse = ", "), ".\n",
      "Got: '", device_type, "'"
    )
  }
  
  if (device_type == "other" && is.null(col_aliases)) {
    stop(
      "load_static_tests(): device_type = 'other' requires `col_aliases` ",
      "mapping all of: date, time, lat, long, height, hdop, sat.\n",
      "Example: col_aliases = list(date = 'Date', time = 'Time', lat = 'Lat', ",
      "long = 'Lon', height = 'Alt_m', hdop = 'HDOP', sat = 'NumSats')"
    )
  }
  
  if (!fs::dir_exists(data_dir)) {
    stop("load_static_tests(): directory not found: ", data_dir)
  }
  
  empty_result <- tibble::tibble(
    static_test_id = character(),
    session_id     = character(),
    device         = character(),
    device_type    = character(),
    field_site     = character(),
    prefix_2       = character(),
    datetime       = as.POSIXct(character()),
    lat            = double(),
    lon            = double(),
    HEIGHT         = double(),
    HDOP           = double()
  )
  
  # -- Step 1 : scan for static test files (both .CSV and .csv) ---------------
  all_files <- unique(c(
    fs::dir_ls(data_dir, glob = "*.CSV", recurse = TRUE),
    fs::dir_ls(data_dir, glob = "*.csv", recurse = TRUE)
  ))
  
  static_files <- all_files[
    stringr::str_detect(
      stringr::str_to_lower(basename(all_files)),
      stringr::str_to_lower(file_pattern)
    )
  ]
  
  if (length(static_files) == 0L) {
    warning(
      "load_static_tests(): no static test files found in '", data_dir, "' ",
      "matching pattern '", file_pattern, "'.\n",
      "Expected filenames like 'staticUM155.CSV' or 'Static_TB03_columbus.csv'.\n",
      "Site-specific thresholds cannot be derived empirically."
    )
    return(empty_result)
  }
  
  cat("load_static_tests(): found", length(static_files),
      "static file(s).\n")
  
  # -- Step 2 : parse filenames -> static_test_id, session_id, prefix_2, device
  #
  # Expected pattern: <anything>[prefix_2][number][optional trailing device label]
  #   "staticUA7"             -> prefix_2 = "UA", static_test_id = "staticUA7",       device = NA
  #   "Static_TN30"           -> prefix_2 = "TN", static_test_id = "Static_TN30",     device = NA
  #   "Static_TB03_columbus"  -> prefix_2 = "TB", static_test_id = "Static_TB03",     device = "columbus"
  #
  # static_test_id is a human-readable id and is NOT guaranteed unique across
  # files -- two device recordings of the same prefix+number share it.
  # session_id IS guaranteed unique per file (it's the full original
  # filename), and is what derive_site_thresholds() uses for per-session
  # centroid grouping, so device-distinct recordings of the same physical
  # test are never pooled together.
  #
  # The 2-letter prefix + following digit run is located via regex; anything
  # after that digit run is treated as a trailing device label (leading "_"
  # or "-" stripped). Matching against prefix_map confirms the prefix is a
  # real site code rather than incidental letters elsewhere in the name.
  file_df <- tibble::tibble(file_path = static_files) |>
    dplyr::mutate(
      fname      = tools::file_path_sans_ext(basename(file_path)),
      session_id = fname,   # full filename -- always unique per file
      
      # Locate [2 letters][1+ digits] anywhere in the name, case-insensitive
      prefix_num_match = stringr::str_extract(
        stringr::str_to_upper(fname), "[A-Z]{2}\\d+"
      ),
      prefix_2_raw = stringr::str_extract(prefix_num_match, "^[A-Z]{2}"),
      
      # Position right after the matched [prefix+digits] block, in the
      # ORIGINAL (not upper-cased) filename, so static_test_id preserves
      # the file's original casing
      match_end = dplyr::if_else(
        is.na(prefix_num_match), NA_integer_,
        stringr::str_locate(stringr::str_to_upper(fname), stringr::fixed(prefix_num_match))[, "end"]
      ),
      
      static_test_id = dplyr::if_else(
        is.na(match_end), NA_character_,
        stringr::str_sub(fname, 1L, match_end)
      ),
      
      device = dplyr::if_else(
        is.na(match_end), NA_character_,
        stringr::str_sub(fname, match_end + 1L) |>
          stringr::str_remove("^[_-]") |>
          (\(x) dplyr::na_if(x, ""))()
      )
    )
  
  # Files where no [2 letters][digits] pattern could be found at all
  no_prefix <- file_df |> dplyr::filter(is.na(prefix_2_raw))
  if (nrow(no_prefix) > 0L) {
    warning(
      "load_static_tests(): cannot extract a site prefix from ", nrow(no_prefix),
      " file(s) -- SKIPPED:\n",
      paste(no_prefix$fname, collapse = "\n")
    )
  }
  
  # Files whose prefix doesn't appear in prefix_map at all
  candidates  <- file_df |> dplyr::filter(!is.na(prefix_2_raw))
  unknown_pfx <- candidates |> dplyr::filter(!prefix_2_raw %in% prefix_map$prefix_2)
  if (nrow(unknown_pfx) > 0L) {
    warning(
      "load_static_tests(): prefix(es) not found in prefix_map -- ",
      nrow(unknown_pfx), " file(s) SKIPPED:\n",
      paste(unknown_pfx$fname, " (prefix '", unknown_pfx$prefix_2_raw, "')",
            sep = "", collapse = "\n"),
      "\nKnown prefixes: ", paste(unique(prefix_map$prefix_2), collapse = ", ")
    )
  }
  
  file_df <- candidates |> dplyr::filter(prefix_2_raw %in% prefix_map$prefix_2)
  
  if (nrow(file_df) == 0L) {
    warning(
      "load_static_tests(): no static files could be matched to a known ",
      "site prefix. Returning an empty result."
    )
    return(empty_result)
  }
  
  # Attach field_site via prefix_map (one row per file -- left_join is safe
  # as long as prefix_map has at most one field_site per prefix_2, which
  # build_prefix_map() already guarantees)
  file_df <- file_df |>
    dplyr::left_join(
      prefix_map |> dplyr::select(prefix_2, field_site),
      by = c("prefix_2_raw" = "prefix_2")
    ) |>
    dplyr::rename(prefix_2 = prefix_2_raw)
  
  # Resolve device type per file (the READER to use -- distinct from the
  # filename-derived `device` label above)
  file_df <- file_df |>
    dplyr::mutate(
      file_device_type = if (device_type == "mix") {
        purrr::map_chr(file_path, .detect_device)
      } else {
        device_type
      }
    )
  
  # -- Step 3 : read and parse each file, dispatching by device type ----------
  .read_one_static <- function(row) {
    
    df <- switch(
      row$file_device_type,
      columbus  = .read_one_columbus(row$file_path),
      catlogger = .read_one_catlogger(row$file_path),
      other     = .read_one_other(row$file_path, col_aliases),
      stop("load_static_tests(): unrecognised device type '", row$file_device_type,
           "' for file '", row$file_path, "'.")
    )
    
    df <- df |>
      dplyr::mutate(
        static_test_id = row$static_test_id,
        session_id     = row$session_id,
        device         = row$device,
        device_type    = row$file_device_type,
        prefix_2       = row$prefix_2,
        field_site     = row$field_site,
        date_str       = stringr::str_pad(as.character(DATE), 6, "left", "0"),
        time_str       = stringr::str_pad(as.character(TIME), 6, "left", "0"),
        datetime       = lubridate::parse_date_time(
          paste(date_str, time_str),
          orders = c("ymdHMS", "ymd HMS"),
          quiet  = TRUE)
      ) |>
      dplyr::select(-date_str, -time_str)
    
    df |>
      dplyr::select(
        static_test_id, session_id, device, device_type, field_site, prefix_2,
        datetime, lat, lon, HEIGHT, HDOP
      )
  }
  
  n_files <- nrow(file_df)
  pb   <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  rows <- purrr::transpose(file_df)
  result_list <- vector("list", n_files)
  
  for (i in seq_len(n_files)) {
    result_list[[i]] <- .read_one_static(rows[[i]])
    utils::setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("\n")
  
  result <- dplyr::bind_rows(result_list)
  
  cat("load_static_tests(): loaded", nrow(result), "GPS fixes from",
      dplyr::n_distinct(result$session_id), "static test session(s) across",
      dplyr::n_distinct(result$static_test_id), "site/number combination(s).\n")
  
  result
}


#' Derive site-specific thresholds from static GPS test data
#'
#' For each field site, computes two thresholds from the static test fixes:
#'
#' - `dist_noise_m`: the 95th percentile of the Haversine distance from the
#'   per-session centroid. Represents the worst-case positional error at
#'   acceptable HDOP values, used as the distance gate in [filter_angle()].
#'
#' - `height_threshold`: the 99.9th percentile of recorded ellipsoidal HEIGHT
#'   plus `height_buffer_m`. Any dog fix above this value cannot represent a
#'   ground-level position and is flagged by [filter_height()].
#'
#' @param static_data A tibble as returned by [load_static_tests()]. Must
#'   contain a `session_id` column (added alongside `static_test_id` so that
#'   device-distinct recordings of the same physical test location are kept
#'   as separate sessions rather than pooled together).
#' @param prefix_map A tibble as returned by [build_prefix_map()], used to
#'   identify which sites have no static test data.
#' @param hdop_threshold Numeric. HDOP ceiling applied to static test fixes' worst bin: 3-5,
#'   before computing `dist_noise_m`. Only fixes with `HDOP <= hdop_threshold`
#'   are used -- this ensures `dist_noise_m` represents the worst-case residual
#'   positional error on fixes that would pass the HDOP filter. Default `4.9`.
#' @param height_overrides A named list of manual `height_threshold` overrides,
#'   keyed by `field_site` (e.g. `list("N'Djamena" = 450)`). Used for sites
#'   without static test files, where you know the correct site-specific
#'   value (e.g. from local elevation data). By default, an entry here is
#'   ONLY applied to a site that has no static test data -- if real static
#'   data already exists for that site, the override is silently ignored
#'   (with a console message noting this) and the static-derived value is
#'   used instead. This protects against an override being filled in
#'   defensively, without knowing whether static data was actually
#'   collected for that site. There is no generic baseline fallback for
#'   height (unlike `dist_noise_m`) -- elevation varies too much from site
#'   to site for a single shared value to be meaningful, even one derived
#'   from an extensive test elsewhere. A site with neither static data nor
#'   an applicable entry here will cause the function to stop (see Details).
#' @param force_height_overrides Logical `TRUE`, or a character vector of
#'   site names, or `character()` (default, meaning none). Explicit opt-in
#'   to let `height_overrides` win even when real static-test data exists
#'   for a site -- use this for the rare case where you know a site's
#'   static test is unrepresentative (e.g. conducted at a different
#'   elevation than the actual deployment area) and want the override to
#'   replace it. Pass `TRUE` to force every site listed in
#'   `height_overrides`, or a character vector to force only specific
#'   sites. Has no effect on sites that have no static data to begin with,
#'   since the override already applies there regardless.
#' @param baseline_dist_noise_m Numeric, or `NULL`. The generic fallback
#'   `dist_noise_m` applied to any site that has no static data. There is no
#'   site-specific override mechanism for `dist_noise_m` -- a site either
#'   has empirical static-test data, or uses this single shared baseline.
#'   This value is not tied to any particular field site -- it's intended
#'   to come from wherever you've collected the most extensive, reliable
#'   static-test data (e.g. a long-running test conducted at a home base
#'   before fieldwork), used as a general-purpose estimate of device noise
#'   for sites where no local static test exists. There is no default;
#'   pass a value explicitly if you want a fallback to be available.
#' @param dog_data Optional tibble as returned by [load_dog_data()]. When
#'   provided, the derived `height_threshold` for each site is validated
#'   against the median dog fix HEIGHT. If the threshold falls below the
#'   median, the function stops with an informative error -- this indicates
#'   that the static tests were conducted at a different elevation than the
#'   dog deployment area. Default `NULL` (no validation).
#'
#' @return A tibble with one row per field site, containing:
#'   - `field_site`       : character
#'   - `dist_noise_m`     : numeric (metres)
#'   - `height_threshold` : numeric (metres, ellipsoidal)
#'   - `dist_source`      : character -- `"static_tests"` or `"baseline"`
#'   - `height_source`    : character -- `"static_tests"` or `"override"`
#'   - `threshold_used`   : character -- a single combined label per row
#'                          summarising whether the baseline fallback was
#'                          used for this site's `dist_noise_m`. One of
#'                          `"static_tests"` (`dist_noise_m` came from
#'                          static data) or `"baseline"` (`dist_noise_m`
#'                          came from `baseline_dist_noise_m`) -- see
#'                          Details. `height_threshold` never contributes to
#'                          this flag, since it has no baseline fallback.
#'
#' @details
#' For sites present in `prefix_map` but absent from `static_data`, the two
#' thresholds are resolved as follows (NOT treated symmetrically):
#'
#' **`height_threshold`** (site-specific overrides only -- no generic
#' baseline, since elevation varies too much site to site for a shared
#' fallback to be meaningful):
#' 1. A site-specific entry in `height_overrides`, if present.
#' 2. Otherwise unresolved (see below) -- there is no fallback step.
#'
#' **`dist_noise_m`** (no site-specific override mechanism -- only real
#' static-test data or the shared baseline):
#' 1. Static-test-derived value, if the site has static data.
#' 2. The generic `baseline_dist_noise_m`, if no static data exists and a
#'    baseline was supplied.
#' 3. Otherwise unresolved (see below).
#'
#' For sites that DO have static data, `height_overrides` is ignored by
#' default -- the static-derived value always wins, even if that site
#' happens to also appear in `height_overrides`. This is deliberate: it
#' protects against an override being filled in defensively, without
#' knowing whether static data exists for that site. Use
#' `force_height_overrides` to explicitly allow the override to win anyway
#' for specific sites where the static test is known to be unrepresentative.
#'
#' Any site left unresolved for either threshold is collected and reported
#' together in a single error at the end of the function, rather than
#' stopping on the first one encountered -- so you can see the full scope of
#' what's missing in one pass.
#'
#' `threshold_used` flags, per site, whether the baseline was invoked for
#' `dist_noise_m` -- this makes it easy to audit afterwards which sites are
#' running on a true empirical distance-noise estimate versus a generic
#' fallback that may not reflect that site's actual device noise
#' characteristics. `height_threshold` is always either static-test-derived
#' or an explicit, deliberately chosen per-site value -- never a guess.
#'
#' Per-session centroids (used to compute `dist_noise_m`) are grouped by
#' `session_id`, not `static_test_id`. This matters when the same site +
#' test number was recorded by more than one device (e.g.
#' `"Static_TB03_columbus"` and `"Static_TB03_catlogger"`): both share
#' `static_test_id = "Static_TB03"` but have distinct `session_id`s, so each
#' device's recording gets its own centroid and is not pooled with the
#' other's positional noise.
#'
#' @examples
#' \dontrun{
#' meta       <- load_metadata("data/metadata.xlsx")
#' pfx_map    <- build_prefix_map(meta)
#' static_df  <- load_static_tests("data/", pfx_map, device_type = "mix")
#'
#' # Site-specific height override for a site with no static data (height
#' # has no generic fallback), plus a distance-noise baseline derived from
#' # an extensive home-base static test, used for any site lacking its own
#' # static data
#' thresholds <- derive_site_thresholds(
#'   static_df, pfx_map,
#'   height_overrides      = list("N'Djamena" = 450),
#'   baseline_dist_noise_m = 35   # e.g. p95 distance noise from a home-base test
#' )
#'
#' # 'Bogo' has static data, but it's known to have been collected at the
#' # wrong elevation -- force the override to replace it for this site only
#' thresholds <- derive_site_thresholds(
#'   static_df, pfx_map,
#'   height_overrides       = list("Bogo" = 380),
#'   force_height_overrides = c("Bogo"),
#'   baseline_dist_noise_m  = 35
#' )
#' print(thresholds)
#' dplyr::filter(thresholds, threshold_used == "baseline")
#' }
#'
#' @export
derive_site_thresholds <- function(static_data,
                                   prefix_map,
                                   dog_data               = NULL,
                                   hdop_threshold         = 4.9,
                                   height_overrides       = list(),
                                   force_height_overrides = character(),
                                   baseline_dist_noise_m  = NULL) {
  
  all_sites <- prefix_map$field_site
  
  # -- Sites with static data -------------------------------------------------
  if (nrow(static_data) > 0L) {
    
    if (!"session_id" %in% names(static_data)) {
      stop(
        "derive_site_thresholds(): `static_data` is missing a `session_id` ",
        "column. This is expected from load_static_tests() and is used to ",
        "group static fixes into per-device, per-test sessions for centroid ",
        "computation. If you constructed `static_data` manually, add a ",
        "`session_id` column (a unique value per recording session)."
      )
    }
    
    empirical <- static_data |>
      dplyr::filter(!is.na(HDOP), HDOP <= hdop_threshold) |>
      dplyr::group_by(field_site, session_id) |>
      dplyr::mutate(
        centroid_lat = mean(lat, na.rm = TRUE),
        centroid_lon = mean(lon, na.rm = TRUE),
        dist_m       = haversine_m(lat, lon, centroid_lat, centroid_lon)
      ) |>
      dplyr::ungroup() |>
      # dist_noise_m = p95 of the worst-case HDOP bin (3--threshold) per site.
      # The worst-case bin drives the distance gate -- using the global p95
      # would underestimate positional error because low-HDOP fixes (small
      # errors, large counts) dominate the distribution and mask the tail.
      dplyr::mutate(
        hdop_bin = cut(HDOP,
                       breaks = c(0, 1, 2, 3, hdop_threshold, Inf),
                       labels = c("<=1", "1-2", "2-3",
                                  paste0("3-", hdop_threshold), "above"))
      ) |>
      dplyr::group_by(field_site, session_id, hdop_bin) |>
      dplyr::summarise(
        p95_dist = quantile(dist_m, 0.95, na.rm = TRUE),
        .groups  = "drop"
      ) |>
      # Take the maximum p95 across all bins and sessions per site
      dplyr::group_by(field_site) |>
      dplyr::summarise(
        dist_noise_m = max(p95_dist, na.rm = TRUE),
        .groups      = "drop"
      )
    
    # Join height_threshold from a separate summarise on the raw static data
    height_summary <- static_data |>
      dplyr::filter(!is.na(HEIGHT)) |>
      dplyr::group_by(field_site) |>
      dplyr::summarise(
        # p99.9 of static HEIGHT = observed vertical scatter ceiling at ground
        # level. No buffer added -- the static fixes ARE the ground truth.
        height_threshold = quantile(HEIGHT, 0.999, na.rm = TRUE),
        .groups = "drop"
      )
    
    empirical <- empirical |>
      dplyr::left_join(height_summary, by = "field_site") |>
      dplyr::mutate(
        dist_noise_m     = round(dist_noise_m),
        height_threshold = round(height_threshold),
        dist_source      = "static_tests",
        height_source    = "static_tests"
      )
    
  } else {
    empirical <- tibble::tibble(
      field_site       = character(),
      dist_noise_m     = double(),
      height_threshold = double(),
      dist_source      = character(),
      height_source    = character()
    )
  }
  
  # -- Guard: validate height_threshold against dog data if provided ----------
  # If dog_data is supplied, check that the derived height_threshold for each
  # site is above the median dog fix HEIGHT. If not, the static tests are not
  # representative of the dog deployment area (different elevation) and the
  # threshold would flag the majority of legitimate fixes.
  
  # Resolve which sites (if any) are allowed to have their height_overrides
  # apply even though real static-test data exists for them. By default,
  # NONE are -- height_overrides only fills gaps where static data is
  # missing, so that an override filled in defensively (without knowing
  # whether static data exists) can never silently shadow real data.
  # Pass force_height_overrides = TRUE to force all sites in height_overrides,
  # or a character vector of specific site names to force just those.
  sites_to_force <- if (isTRUE(force_height_overrides)) {
    names(height_overrides)
  } else if (is.character(force_height_overrides)) {
    force_height_overrides
  } else {
    character()
  }
  
  if (!is.null(dog_data)) {
    
    dog_height_summary <- dog_data |>
      dplyr::group_by(field_site) |>
      dplyr::summarise(
        median_dog_height = median(HEIGHT, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Apply height_overrides BEFORE checking, but ONLY for sites in
    # sites_to_force -- a site whose static-derived threshold is too low
    # but has a FORCED override must not trigger the error. A site with an
    # override that was NOT forced should still be checked against its real
    # static-derived value, since that override will be ignored later anyway.
    empirical_to_check <- empirical |>
      dplyr::mutate(
        height_threshold = dplyr::if_else(
          field_site %in% sites_to_force,
          purrr::map_dbl(field_site, ~ {
            if (.x %in% sites_to_force && .x %in% names(height_overrides))
              as.double(height_overrides[[.x]])
            else
              NA_real_
          }),
          height_threshold
        )
      )
    
    check <- empirical_to_check |>
      dplyr::left_join(dog_height_summary, by = "field_site") |>
      dplyr::filter(!is.na(median_dog_height)) |>
      dplyr::filter(height_threshold < median_dog_height,
                    !field_site %in% sites_to_force)
    
    if (nrow(check) > 0L) {
      problem_lines <- purrr::map_chr(seq_len(nrow(check)), function(i) {
        sprintf(
          "  '%s': height_threshold = %d m < median dog HEIGHT = %d m",
          check$field_site[i],
          as.integer(check$height_threshold[i]),
          as.integer(check$median_dog_height[i])
        )
      })
      stop(
        "derive_site_thresholds(): the height_threshold derived from static ",
        "tests is BELOW the median dog fix HEIGHT for the following site(s):\n",
        paste(problem_lines, collapse = "\n"), "\n\n",
        "This means the static tests were conducted at a different elevation ",
        "than the dog deployment area and cannot be used to derive a valid ",
        "height_threshold.\n\n",
        "Solution: provide a manual override AND force it (since static data ",
        "already exists for this site, the override is ignored unless ",
        "forced):\n",
        paste(
          sprintf("  height_overrides = list('%s' = <value>)\n  force_height_overrides = c('%s')",
                  check$field_site, check$field_site),
          collapse = "\n"
        ), "\n",
        "Tip: use the median dog HEIGHT + 250 m as a starting point:\n",
        paste(
          sprintf("  '%s' = %d",
                  check$field_site,
                  as.integer(check$median_dog_height + 250)),
          collapse = "\n"
        )
      )
    }
  }
  
  sites_with_data <- empirical$field_site
  sites_missing   <- setdiff(all_sites, sites_with_data)
  
  # -- Sites without static data -- resolve via override (height) / baseline
  #    (distance), then collect as unresolved if neither applies -----------
  # height_threshold: site-specific override -> unresolved (no baseline)
  # dist_noise_m:     (no override mechanism) -> baseline -> unresolved
  unresolved_sites <- character()
  
  override_rows <- purrr::map_dfr(sites_missing, function(site) {
    
    has_height_override <- site %in% names(height_overrides)
    
    height_resolved <- has_height_override
    dist_resolved    <- !is.null(baseline_dist_noise_m)
    
    if (!height_resolved || !dist_resolved) {
      unresolved_sites <<- c(unresolved_sites, site)
      return(NULL)   # skip row construction; this site will error out below
    }
    
    height_value  <- as.double(height_overrides[[site]])
    height_source <- "override"
    
    dist_value  <- as.double(baseline_dist_noise_m)
    dist_source <- "baseline"
    
    cat("derive_site_thresholds(): site '", site, "' has no static tests; ",
        "height_threshold from override, dist_noise_m from baseline.\n",
        sep = "")
    
    tibble::tibble(
      field_site       = site,
      dist_noise_m     = dist_value,
      height_threshold = height_value,
      dist_source      = dist_source,
      height_source    = height_source
    )
  })
  
  if (length(unresolved_sites) > 0L) {
    problem_lines <- purrr::map_chr(unresolved_sites, function(site) {
      missing_pieces <- c(
        if (!(site %in% names(height_overrides)))
          "height_threshold (no override -- height has no baseline fallback)",
        if (is.null(baseline_dist_noise_m))
          "dist_noise_m (no baseline_dist_noise_m)"
      )
      sprintf("  '%s': missing %s", site, paste(missing_pieces, collapse = " and "))
    })
    stop(
      "derive_site_thresholds(): ", length(unresolved_sites),
      " site(s) have no static tests and cannot be resolved:\n",
      paste(problem_lines, collapse = "\n"), "\n\n",
      "For height_threshold, provide a site-specific override ",
      "(height_overrides) -- there is no generic baseline for height, since ",
      "elevation varies too much from site to site.\n",
      "For dist_noise_m, there is no per-site override -- provide ",
      "baseline_dist_noise_m to cover any site without static data.\n",
      "Tip: local elevation (m a.s.l.) + ~150 m is a reasonable starting ",
      "point for height_threshold."
    )
  }
  
  result <- dplyr::bind_rows(empirical, override_rows) |>
    dplyr::arrange(field_site)
  
  # Apply height_overrides to sites that ALREADY have a static-test-derived
  # height_threshold (i.e. sites in `empirical`, not the override_rows we
  # just built for missing sites -- those already used the override value).
  # By default this does NOT happen: if real static data exists for a site,
  # it wins, even if that site also happens to appear in height_overrides --
  # this protects against someone filling in an override defensively without
  # knowing whether static data was actually collected for that site.
  # sites_to_force (computed earlier) is the explicit opt-in for the rare
  # case where the static test IS known to be unrepresentative (e.g. wrong
  # elevation) and the override should win anyway.
  ignored_overrides <- setdiff(
    intersect(names(height_overrides), sites_with_data),
    sites_to_force
  )
  if (length(ignored_overrides) > 0L) {
    cat("derive_site_thresholds(): height_overrides for the following site(s) ",
        "were IGNORED because static-test data already exists for them: ",
        paste(ignored_overrides, collapse = ", "), ".\n",
        "If you intended to override real static data (e.g. because the ",
        "test is known to be unrepresentative), pass these site names via ",
        "`force_height_overrides`.\n", sep = "")
  }
  
  if (length(sites_to_force) > 0L) {
    result <- result |>
      dplyr::mutate(
        height_threshold = purrr::map_dbl(field_site, ~ {
          if (.x %in% sites_to_force && .x %in% names(height_overrides))
            as.double(height_overrides[[.x]])
          else
            height_threshold[field_site == .x]
        }),
        height_source = dplyr::if_else(
          field_site %in% sites_to_force & field_site %in% names(height_overrides),
          "override",
          height_source
        )
      )
  }
  
  # threshold_used: per-site summary flag -- "baseline" if dist_noise_m for
  # that site came from the generic baseline, "static_tests" otherwise.
  # height_source can only ever be "static_tests" or "override" (height has
  # no baseline fallback), so it never contributes to this flag.
  result <- result |>
    dplyr::mutate(
      threshold_used = dplyr::if_else(
        dist_source == "baseline",
        "baseline",
        "static_tests"
      )
    )
  
  # -- Print summary ----------------------------------------------------------
  cat("\n--- Site thresholds ---\n")
  print(result, n = Inf)
  cat("\n")
  
  result
}