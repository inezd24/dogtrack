# =============================================================================
# R/load_data.R
# Package : dogtrack
# Purpose : Scan a directory for dog GPS CSV files from one or more supported
#           device types, apply session-A selection, parse and normalise
#           each file, then join with metadata.
#
# Exported functions:
#   load_dog_data() -- returns a single tidy tibble, one row per GPS fix,
#                     with metadata columns joined and field_site attached
#
# Internal helpers (not exported):
#   .detect_device()      -- per-file auto-detection used when device_type = "mix"
#   .read_one_columbus()   -- Columbus P-10 Pro reader
#   .read_one_catlogger()  -- CatLog reader
#   .read_one_other()      -- generic reader driven by user-supplied col_aliases
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Detects whether a CSV file was produced by a Columbus P-10 Pro or a CatLog
# device, based on structural signatures rather than the file extension:
#   - CatLog files have the literal text "Name:CatLog" in the first cell
#   - Columbus files have INDEX and TAG columns in row 1
# Only used when device_type = "mix". Returns "columbus" or "catlogger"; stops
# with an informative error if neither signature is found (since "other"
# devices cannot be auto-detected -- they must be called separately).
# Assumes the file has already been checked for emptiness by the caller --
# this function does not itself guard against 0-row files, beyond the
# explicit check below which converts that case into a clear error message.
# -----------------------------------------------------------------------------
.detect_device <- function(file_path) {
  
  first_cell <- suppressMessages(
    readr::read_csv(file_path, col_names = FALSE, n_max = 1,
                    col_types = readr::cols(.default = "c"),
                    show_col_types = FALSE)
  )
  
  if (nrow(first_cell) == 0L || ncol(first_cell) == 0L) {
    stop(
      "device auto-detection failed: file '", basename(file_path), "' ",
      "appears to be empty and cannot be used for device auto-detection."
    )
  }
  
  first_cell_value <- first_cell[[1, 1]]
  
  if (!is.na(first_cell_value) &&
      stringr::str_detect(first_cell_value, stringr::regex("Name:CatLog", ignore_case = TRUE))) {
    return("catlogger")
  }
  
  header_row1 <- suppressMessages(
    readr::read_csv(file_path, n_max = 0, show_col_types = FALSE)
  )
  if (all(c("INDEX", "TAG") %in% names(header_row1))) {
    return("columbus")
  }
  
  stop(
    "device auto-detection failed for '", basename(file_path), "' ",
    "(device_type = 'mix').\n",
    "Expected either a CatLog file (\"Name:CatLog\" in the first cell) or a ",
    "Columbus file (INDEX and TAG columns in row 1).\n",
    "If this file is from a different device, load it separately with ",
    "device_type = 'other' and col_aliases."
  )
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Reads a single Columbus P-10 Pro CSV. Header is in row 1. LATITUDE/LONGITUDE
# arrive as "LATITUDE N/S" / "LONGITUDE E/W" and need sign conversion
# (handled inside preprocess_gps(), unchanged from the original pipeline).
# INDEX is left untyped (not forced via col_types) since readr::col_guess()
# is not a valid named-column type and triggers a spurious parsing warning;
# INDEX is not used downstream so its inferred type doesn't matter.
# -----------------------------------------------------------------------------
.read_one_columbus <- function(file_path) {
  
  df <- readr::read_csv(
    file_path,
    col_types = readr::cols(
      TAG             = readr::col_character(),
      DATE            = readr::col_character(),
      TIME            = readr::col_character(),
      `LATITUDE N/S`  = readr::col_character(),
      `LONGITUDE E/W` = readr::col_character(),
      HEIGHT          = readr::col_double(),
      SPEED           = readr::col_double(),
      HEADING         = readr::col_double(),
      SAT             = readr::col_integer(),
      HDOP            = readr::col_double(),
      .default        = readr::col_guess()   # INDEX and any other extras inferred normally
    ),
    show_col_types = FALSE
  ) |>
    preprocess_gps()  # existing pipeline step: resolves lat/lon sign, etc.
  
  df
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Reads a single CatLog CSV. The file has a device banner in row 1
# ("Name:CatLog") with the real header starting on row 7. Columns arrive
# capitalised but not in SCREAMING_SNAKE like Columbus, so they're
# snake_cased and renamed to the canonical schema directly (no col_aliases
# needed -- CatLog's layout is fixed and known).
# CatLog has no INDEX column (one is generated here as the row number) and
# no HEADING column (set to NA, since the device doesn't report it).
# Lat/long arrive as plain decimal degrees, so no N/S, E/W sign conversion
# is needed (unlike Columbus).
# -----------------------------------------------------------------------------
.read_one_catlogger <- function(file_path) {
  
  df <- suppressMessages(
    readr::read_csv(
      file_path,
      skip = 6,   # banner + blank rows + "--------" occupy rows 1-6; header is row 7
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )
  )
  
  names(df) <- stringr::str_remove(
    stringr::str_replace_all(
      stringr::str_replace_all(stringr::str_to_lower(names(df)), "[^a-z0-9]", "_"),
      "_+", "_"),
    "^_|_$"
  )
  
  required_catlogger <- c(
    date = "date", time = "time", lat = "latitude", long = "longitude",
    height = "altitude", hdop = "hdop", sat = "satellites"
  )
  missing <- setdiff(required_catlogger, names(df))
  if (length(missing) > 0) {
    stop(
      "load_dog_data(): CatLog file '", basename(file_path), "' is missing ",
      "expected column(s): ", paste(missing, collapse = ", "), ".\n",
      "Columns found: ", paste(names(df), collapse = ", ")
    )
  }
  
  df <- df |>
    dplyr::rename(
      DATE   = date,
      TIME   = time,
      HEIGHT = altitude,
      HDOP   = hdop,
      SAT    = satellites
    ) |>
    dplyr::mutate(
      INDEX   = dplyr::row_number(),
      TAG     = NA_character_,
      HEIGHT  = as.numeric(HEIGHT),
      SAT     = as.integer(SAT),
      HDOP    = as.numeric(HDOP),
      SPEED   = suppressWarnings(as.numeric(`speed_km_h`)),
      HEADING = NA_real_,
      lat     = suppressWarnings(as.numeric(latitude)),
      lon     = suppressWarnings(as.numeric(longitude))
    )
  
  df
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Reads a single CSV from an unsupported ("other") device, using a
# user-supplied col_aliases mapping to locate the seven required columns:
# date, time, lat, long, height, hdop, sat. Unlike load_metadata()'s
# col_aliases, ALL seven must be supplied explicitly every time -- there is
# no built-in variant table for arbitrary devices, since their column
# conventions are unknown to the package.
# Assumes header is in row 1 with no banner rows; if a particular device has
# banner rows above its header, the user should pre-process the file before
# calling load_dog_data(), since there is no generic way to detect where an
# unknown device's header begins.
# -----------------------------------------------------------------------------
.read_one_other <- function(file_path, col_aliases) {
  
  required_other <- c("date", "time", "lat", "long", "height", "hdop", "sat")
  
  if (is.null(col_aliases) || !is.list(col_aliases) || is.null(names(col_aliases))) {
    stop(
      "load_dog_data(): device_type = 'other' requires `col_aliases`, a named ",
      "list mapping ALL of the following required columns to your file's ",
      "actual column names:\n",
      "  ", paste(required_other, collapse = ", "), "\n",
      "Example:\n",
      "  col_aliases = list(date = 'Date', time = 'Time', lat = 'Lat',\n",
      "                      long = 'Lon', height = 'Alt_m', hdop = 'HDOP',\n",
      "                      sat = 'NumSats')"
    )
  }
  
  missing_keys <- setdiff(required_other, names(col_aliases))
  if (length(missing_keys) > 0) {
    stop(
      "load_dog_data(): `col_aliases` is missing required key(s): ",
      paste(missing_keys, collapse = ", "), ".\n",
      "device_type = 'other' requires ALL of: ", paste(required_other, collapse = ", "),
      " to be mapped explicitly."
    )
  }
  
  df_raw <- suppressMessages(
    readr::read_csv(file_path, col_types = readr::cols(.default = "c"),
                    show_col_types = FALSE)
  )
  
  missing_cols <- setdiff(unlist(col_aliases[required_other]), names(df_raw))
  if (length(missing_cols) > 0) {
    stop(
      "load_dog_data(): file '", basename(file_path), "' is missing column(s) ",
      "named in `col_aliases`: ", paste(missing_cols, collapse = ", "), ".\n",
      "Columns found in file: ", paste(names(df_raw), collapse = ", ")
    )
  }
  
  df <- df_raw |>
    dplyr::rename(
      DATE = !!col_aliases[["date"]],
      TIME = !!col_aliases[["time"]],
      lat  = !!col_aliases[["lat"]],
      lon  = !!col_aliases[["long"]],
      HEIGHT = !!col_aliases[["height"]],
      HDOP   = !!col_aliases[["hdop"]],
      SAT    = !!col_aliases[["sat"]]
    ) |>
    dplyr::mutate(
      INDEX   = dplyr::row_number(),
      TAG     = NA_character_,
      SPEED   = NA_real_,
      HEADING = NA_real_,
      lat     = suppressWarnings(as.numeric(lat)),
      lon     = suppressWarnings(as.numeric(lon)),
      HEIGHT  = suppressWarnings(as.numeric(HEIGHT)),
      HDOP    = suppressWarnings(as.numeric(HDOP)),
      SAT     = suppressWarnings(as.integer(SAT))
    )
  
  df
}


#' Load and merge dog GPS CSV files with metadata
#'
#' Scans a directory recursively for dog GPS CSV files, excluding files whose
#' name begins with `"static"`. Supports Columbus P-10 Pro and CatLog devices
#' natively, plus any other device via an explicit column mapping. Applies
#' session-A selection logic to handle dogs with multiple recording sessions,
#' detects file-naming conflicts caused by appended device numbers, parses
#' datetime and coordinates, computes `delta_t_sec`, and joins the result
#' with the project metadata.
#'
#' @param data_dir Character. Path to the directory containing the dog GPS
#'   CSV files (static test files in the same folder are automatically
#'   excluded, as are empty files and files that don't match the expected
#'   naming convention, with a warning in each case).
#' @param metadata A tibble as returned by [load_metadata()].
#' @param device_type Character. One of:
#'   \describe{
#'     \item{`"columbus"`}{All files in `data_dir` are read as Columbus
#'       P-10 Pro output.}
#'     \item{`"catlogger"`}{All files in `data_dir` are read as CatLog output.}
#'     \item{`"other"`}{All files in `data_dir` are read using the column
#'       mapping supplied in `col_aliases`, which is required in this case.}
#'     \item{`"mix"`}{Each file's device type is detected automatically
#'       (Columbus vs CatLog only -- `"other"`-type files cannot be mixed
#'       into auto-detection and must be loaded in a separate call).}
#'   }
#' @param col_aliases Named list, or `NULL`. Required when `device_type =
#'   "other"`. Must map ALL of the following canonical names to the column
#'   names actually used in your files: `date`, `time`, `lat`, `long`,
#'   `height`, `hdop`, `sat`. Ignored for `"columbus"`, `"catlogger"`, and
#'   `"mix"`, since those use fixed, known layouts.
#'
#' @return A tibble with one row per GPS fix, containing all original device
#'   columns plus:
#'   - `file_key`     : character -- full filename without extension, unique
#'                      per CSV file (e.g. `"UMR001-01A"`)
#'   - `join_key`     : character -- canonical base key used for metadata
#'                      join (e.g. `"TUA011-01"`), with any session-letter or
#'                      device-number/commentary suffix stripped
#'   - `version`      : character, or `NA` -- the device-number suffix
#'                      (e.g. `"28"`), populated ONLY when multiple files
#'                      share the same `join_key`; otherwise `NA`. Lets you
#'                      distinguish and filter between conflicting files for
#'                      the same dog deployment without altering `join_key`.
#'   - `id_comment`   : character, or `NA` -- any free-text commentary found
#'                      after the device number (e.g. `"not fully sure about
#'                      the ID"`), populated under the same condition as
#'                      `version`.
#'   - `lat`          : numeric -- decimal degrees (S = negative)
#'   - `lon`          : numeric -- decimal degrees (W = negative)
#'   - `datetime`     : POSIXct -- parsed from DATE + TIME columns
#'   - `delta_t_sec`  : numeric -- seconds elapsed since previous fix
#'                      within the same dog track; `NA` for the first fix
#'   - `prefix_2`     : character -- 2-letter site prefix (e.g. `"UM"`)
#'   - `setting`      : character -- environment code extracted from filename
#'                      (e.g. `"R"` for rural, `"U"` for urban); `NA` if
#'                      absent
#'   - `device_type`  : character -- `"columbus"`, `"catlogger"`, or `"other"`,
#'                      recorded per file for traceability
#'   - All metadata columns joined from [load_metadata()]
#'
#' @details
#' **Device support:** Columbus and CatLog devices are read with built-in,
#' fixed column layouts. Columbus files have a single header row; CatLog
#' files have a banner ("Name:CatLog") with the real header starting on row 7.
#' CatLog files have no native `INDEX` or `TAG` column (a row-number index is
#' generated, and `TAG` is set to `NA`) and no `HEADING` column (set to `NA`).
#' CatLog's latitude/longitude arrive as plain decimal degrees and need no
#' sign conversion, unlike Columbus's `"LATITUDE N/S"` / `"LONGITUDE E/W"`
#' format.
#'
#' For any other device, set `device_type = "other"` and supply `col_aliases`
#' mapping all seven required columns (`date`, `time`, `lat`, `long`,
#' `height`, `hdop`, `sat`) to your file's actual column names. Unlike
#' [load_metadata()]'s `col_aliases`, there is no built-in variant table here
#' -- every required column must be mapped explicitly, and the function stops
#' with an error if any are missing. `"other"` devices are assumed to have a
#' header in row 1 with no banner rows above it; if your device has banner
#' rows, pre-process the file before calling this function.
#'
#' **Mixed batches:** set `device_type = "mix"` to auto-detect Columbus vs
#' CatLog per file (CatLog is identified by `"Name:CatLog"` in the first
#' cell; Columbus by `INDEX`/`TAG` columns in row 1). `"other"`-type files
#' cannot be included in a mixed batch -- call the function separately for
#' those with `device_type = "other"`.
#'
#' **File discovery is case-insensitive for the extension:** both `.CSV` and
#' `.csv` files are matched, since naming conventions vary across sites and
#' devices, and some filesystems treat the extension case-sensitively.
#'
#' **Filename convention and parsing:** filenames are expected to follow
#' `[Country][Location][Setting]-[house: 3 digits]-[individual: 2 digits]`
#' (e.g. `"TUA011-01"`). The separators may be `-` or `_` (e.g.
#' `"UAR001_01.CSV"` is accepted and normalised to `"UAR001-01"` for
#' `join_key`), since both conventions appear across sites/eras. This is
#' optionally followed by one of two mutually exclusive trailing suffixes:
#' \itemize{
#'   \item a session letter (`A`/`B`/`C`), glued directly or with one
#'     leading hyphen (e.g. `"TUA011-01A"` or `"TUA011-01-A"`)
#'   \item a device number, introduced by a hyphen, optionally followed by
#'     further hyphen-or-space-separated free text (e.g. `"TUA011-01-28"` or
#'     `"TUA011-01-102-not fully sure about the ID"`)
#' }
#' The canonical core (without either suffix) becomes `join_key`. The device
#' number and any trailing commentary become `version` and `id_comment`
#' respectively, but ONLY when more than one file shares the same canonical
#' core -- a single file with a device-number suffix and no sibling files
#' has `version` and `id_comment` set to `NA`, since there's nothing to
#' disambiguate. This ensures `join_key` (and therefore the dog's identity)
#' is never altered by these suffixes, regardless of how many files exist
#' for it.
#'
#' **Files that don't match the naming convention at all** (e.g. a name with
#' no parseable house/individual number, such as `"TNR-76-or-79.csv"`) cannot
#' be assigned a `join_key` and are excluded automatically, with a warning
#' listing the affected filenames. These typically indicate genuine
#' ambiguity in the field record (e.g. uncertainty about which of two dogs a
#' file belongs to) and need to be resolved manually -- by renaming the file
#' once the correct identity is confirmed -- before they can be loaded.
#'
#' **Session selection:** some dogs have multiple CSV files for the same
#' deployment (suffixes `A`, `B`, `C` or `1`, `2` in the filename). When a
#' file with suffix `A` exists for a given base key, all other sessions
#' for that key are dropped. When no `A` suffix exists, the single
#' unsuffixed file is kept. This logic is independent of `version` --
#' files distinguished only by device number (no session letters at all)
#' are NOT subject to session-A filtering and are all kept.
#'
#' **Empty files:** any CSV file that is 0 bytes is skipped automatically,
#' with a warning listing the skipped file paths, before any reading or
#' device-detection is attempted.
#'
#' **Unmatched files:** CSV files with no metadata match are excluded via
#' inner join and reported as a warning. This is a data quality flag --
#' investigate before proceeding.
#'
#' **Datetime parsing:** DATE and TIME may arrive as human-readable strings
#' (`"2026-01-07"`, `"09:47:09"`) or compact numerics (`260107`, `94709`).
#' Both are handled automatically.
#'
#' **Progress:** a text progress bar (`utils::txtProgressBar`) is shown while
#' files are read and parsed, since this can be slow for large batches.
#'
#' @examples
#' \dontrun{
#' meta <- load_metadata("data/metadata.xlsx")
#'
#' # All files in the directory are Columbus output (original behaviour)
#' dog_data <- load_dog_data("data/", meta, device_type = "columbus")
#'
#' # All files are CatLog output
#' dog_data <- load_dog_data("data/", meta, device_type = "catlogger")
#'
#' # Mixed folder of Columbus and CatLog files
#' dog_data <- load_dog_data("data/", meta, device_type = "mix")
#'
#' # An unsupported device
#' dog_data <- load_dog_data(
#'   "data/", meta, device_type = "other",
#'   col_aliases = list(date = "Date", time = "Time", lat = "Lat",
#'                       long = "Lon", height = "Alt_m", hdop = "HDOP",
#'                       sat = "NumSats")
#' )
#'
#' # Inspect files that needed disambiguation
#' dog_data |> dplyr::filter(!is.na(version)) |> dplyr::distinct(join_key, version, id_comment)
#' }
#'
#' @export
load_dog_data <- function(data_dir, metadata, device_type, col_aliases = NULL) {
  
  valid_device_types <- c("columbus", "catlogger", "other", "mix")
  if (!device_type %in% valid_device_types) {
    stop(
      "load_dog_data(): `device_type` must be one of: ",
      paste(valid_device_types, collapse = ", "), ".\n",
      "Got: '", device_type, "'"
    )
  }
  
  if (device_type == "other" && is.null(col_aliases)) {
    stop(
      "load_dog_data(): device_type = 'other' requires `col_aliases` mapping ",
      "all of: date, time, lat, long, height, hdop, sat.\n",
      "Example: col_aliases = list(date = 'Date', time = 'Time', lat = 'Lat', ",
      "long = 'Lon', height = 'Alt_m', hdop = 'HDOP', sat = 'NumSats')"
    )
  }
  
  if (!fs::dir_exists(data_dir)) {
    stop("load_dog_data(): directory not found: ", data_dir)
  }
  
  # -- Step 1 : scan for dog CSV files (exclude static tests) -----------------
  # Match both .CSV and .csv -- glob matching is case-sensitive on some
  # filesystems, and naming conventions vary across sites/devices, so both
  # extensions are matched explicitly rather than relying on glob alone.
  all_files <- unique(c(
    fs::dir_ls(data_dir, glob = "*.CSV", recurse = TRUE),
    fs::dir_ls(data_dir, glob = "*.csv", recurse = TRUE)
  ))
  
  dog_files <- all_files[
    !stringr::str_detect(
      stringr::str_to_lower(basename(all_files)), "^static"
    )
  ]
  
  if (length(dog_files) == 0L) {
    stop(
      "load_dog_data(): no dog GPS CSV files found in '", data_dir, "'.\n",
      "Check that the path is correct and files have a .csv/.CSV extension."
    )
  }
  
  # -- Step 2 : parse filenames -> file_key, base_key, session, version ------
  
  # Drop empty files before any further processing, with a warning per file.
  # An empty CSV has no header and no data to detect a device from or parse,
  # regardless of device_type.
  file_sizes  <- fs::file_size(dog_files)
  empty_files <- dog_files[file_sizes == 0]
  
  if (length(empty_files) > 0L) {
    warning(
      "load_dog_data(): ", length(empty_files),
      " file(s) are empty and will be SKIPPED:\n",
      paste(empty_files, collapse = "\n")
    )
    dog_files <- dog_files[file_sizes > 0]
  }
  
  if (length(dog_files) == 0L) {
    stop(
      "load_dog_data(): no non-empty dog GPS CSV files found in '", data_dir, "'.\n",
      "Check that the path is correct and files have a .csv/.CSV extension."
    )
  }
  
  # Filename convention: [Country][Location][Setting]-[house:3 digits]-[individual:2 digits]
  #   e.g. "TUA011-01" is the canonical core.
  # The separator between the 3-letter site code and the house number, and
  # between the house number and individual number, may be "-" OR "_"
  # (older Uganda files use "_", e.g. "UAR001_01.CSV"); both are normalised
  # to "-" in `base_key` so join_key is always consistent regardless of the
  # original file's separator style.
  # Two kinds of trailing suffix are recognised, and they are mutually exclusive:
  #   (a) a session letter (A/B/C), glued with or without a hyphen/underscore
  #       e.g. "TUA011-01A" or "TUA011-01-A" -> session = "A"
  #   (b) a device number, optionally followed by free-text commentary,
  #       always introduced by a hyphen after the canonical core
  #       e.g. "TUA011-01-28" -> device_suffix = "28"
  #            "TUA011-01-102-not fully sure about the ID"
  #              -> device_suffix = "102", id_comment_raw = "not fully sure about the ID"
  # `version` and `id_comment` are only populated downstream when multiple
  # files share the same base_key -- a lone file's device number/commentary
  # is discarded, since there's nothing to disambiguate.
  #
  # Files with no parseable core (e.g. "TNR-76-or-79.csv", where the field
  # team recorded genuine uncertainty about device identity directly in the
  # filename) cannot be assigned a join_key and are excluded separately below.
  file_df <- tibble::tibble(file_path = dog_files) |>
    dplyr::mutate(
      file_key = tools::file_path_sans_ext(basename(file_path)),
      
      # Canonical core: COUNTRY+LOCATION+SETTING(3 letters) [-_]? house(3 digits) [-_] individual(2 digits)
      # Captured as a single match, then normalised to "-" separators below.
      core_raw = stringr::str_extract(file_key, "^[A-Z]{3}[-_]?\\d{3}[-_]\\d{2}"),
      
      # Normalise to a canonical "-"-separated form regardless of which
      # separator the original filename used, so base_key/join_key is
      # always consistent (e.g. "UAR001_01" and "UAR001-01" both become
      # the same base_key if they ever both occurred).
      core = dplyr::if_else(
        is.na(core_raw), NA_character_,
        stringr::str_replace_all(core_raw, "_", "-")
      ),
      
      # Everything after the ORIGINAL (unnormalised) core match, since
      # remainder must be sliced from the real filename, not the
      # normalised version
      remainder = dplyr::if_else(
        is.na(core_raw), NA_character_,
        stringr::str_sub(file_key, stringr::str_length(core_raw) + 1L)
      ),
      
      # (a) Session letter: glued directly, or with one leading hyphen, then nothing else
      session = dplyr::case_when(
        is.na(remainder) ~ NA_character_,
        stringr::str_detect(remainder, "^[A-Za-z]$")  ~ stringr::str_to_upper(remainder),
        stringr::str_detect(remainder, "^-[A-Za-z]$") ~ stringr::str_to_upper(stringr::str_sub(remainder, 2, 2)),
        TRUE ~ NA_character_
      ),
      
      # (b) Device number + optional trailing commentary, only when (a) didn't match
      #     and remainder starts with "-<digits>"
      device_suffix = dplyr::if_else(
        !is.na(remainder) & is.na(session) & stringr::str_detect(remainder, "^-\\d+"),
        stringr::str_extract(remainder, "^-\\d+") |> stringr::str_remove("^-"),
        NA_character_
      ),
      id_comment_raw = dplyr::if_else(
        !is.na(device_suffix),
        remainder |>
          stringr::str_remove("^-\\d+") |>   # drop the device number chunk
          stringr::str_remove("^-") |>       # drop the separating hyphen, if present
          (\(x) dplyr::na_if(x, ""))(),       # empty string -> NA
        NA_character_
      ),
      
      base_key = core,   # the canonical core IS the base_key now -- nothing else to strip
      
      prefix_2 = stringr::str_extract(basename(file_path), "^[A-Z]{2}"),
      setting  = stringr::str_extract(basename(file_path), "^[A-Z]{3}") |>
        stringr::str_sub(3L, 3L)
    )
  
  # Files that don't match the naming convention at all (core = NA) cannot be
  # assigned a join_key. Exclude them with a clear warning rather than letting
  # them silently collapse into a shared NA group.
  unparseable <- file_df |> dplyr::filter(is.na(core))
  if (nrow(unparseable) > 0L) {
    warning(
      "load_dog_data(): ", nrow(unparseable),
      " file(s) do not match the expected naming convention ",
      "([Country][Location][Setting]-[house]-[individual]) and are EXCLUDED:\n",
      paste(basename(unparseable$file_path), collapse = "\n"),
      "\nThese filenames cannot be reliably assigned to a dog. ",
      "Resolve the correct identity and rename the file before re-running."
    )
    file_df <- file_df |> dplyr::filter(!is.na(core))
  }
  
  if (nrow(file_df) == 0L) {
    stop(
      "load_dog_data(): no files in '", data_dir, "' match the expected ",
      "naming convention after exclusions. Nothing to load."
    )
  }
  
  # Resolve device type per file (only AFTER unparseable files are dropped,
  # since there's no point running detection on files we'll exclude anyway)
  file_df <- file_df |>
    dplyr::mutate(
      file_device_type = if (device_type == "mix") {
        purrr::map_chr(file_path, .detect_device)
      } else {
        device_type
      }
    ) |>
    dplyr::group_by(base_key) |>
    dplyr::mutate(
      has_session_A   = any(session == "A", na.rm = TRUE),
      n_with_base_key = dplyr::n(),
      
      keep = dplyr::case_when(
        session == "A"                  ~ TRUE,
        is.na(session) & !has_session_A ~ TRUE,
        TRUE                            ~ FALSE
      ),
      
      version    = dplyr::if_else(n_with_base_key > 1L, device_suffix,  NA_character_),
      id_comment = dplyr::if_else(n_with_base_key > 1L, id_comment_raw, NA_character_)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(keep) |>
    dplyr::select(-core, -core_raw, -remainder, -device_suffix, -id_comment_raw,
                  -has_session_A, -n_with_base_key, -keep)
  
  cat("load_dog_data(): found", nrow(file_df), "dog CSV file(s) after session selection.\n")
  if (device_type == "mix") {
    cat("load_dog_data(): device types detected --",
        paste(names(table(file_df$file_device_type)), table(file_df$file_device_type),
              sep = ": ", collapse = ", "), "\n")
  }
  if (any(!is.na(file_df$version))) {
    n_versioned <- sum(!is.na(file_df$version))
    cat("load_dog_data():", n_versioned,
        "file(s) disambiguated via `version` (multiple files shared a base_key).\n")
  }
  
  # -- Step 3 : read and parse each CSV, dispatching by device type -----------
  .read_one <- function(row) {
    
    df <- switch(
      row$file_device_type,
      columbus  = .read_one_columbus(row$file_path),
      catlogger = .read_one_catlogger(row$file_path),
      other     = .read_one_other(row$file_path, col_aliases),
      stop("load_dog_data(): unrecognised device type '", row$file_device_type,
           "' for file '", row$file_path, "'.")
    )
    
    df <- df |>
      dplyr::mutate(
        file_key    = row$file_key,
        join_key    = row$base_key,
        version     = row$version,
        id_comment  = row$id_comment,
        prefix_2    = row$prefix_2,
        setting     = row$setting,
        device_type = row$file_device_type,
        date_str    = stringr::str_pad(as.character(DATE), 6, "left", "0"),
        time_str    = stringr::str_pad(as.character(TIME), 6, "left", "0"),
        datetime    = lubridate::parse_date_time(
          paste(date_str, time_str),
          orders = c("ymdHMS", "ymd HMS"),
          quiet  = TRUE)
      ) |>
      dplyr::select(-date_str, -time_str) |>
      dplyr::arrange(datetime) |>
      dplyr::mutate(
        delta_t_sec = as.numeric(
          difftime(datetime, dplyr::lag(datetime), units = "secs")
        )
      )
    
    df
  }
  
  n_files <- nrow(file_df)
  pb   <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  rows <- purrr::transpose(file_df)
  dog_data_list <- vector("list", n_files)
  
  for (i in seq_len(n_files)) {
    dog_data_list[[i]] <- .read_one(rows[[i]])
    utils::setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("\n")  # progress bar leaves the cursor on the same line; move past it
  
  dog_data_raw <- dplyr::bind_rows(dog_data_list)
  
  cat("load_dog_data(): loaded", nrow(dog_data_raw), "GPS fixes from",
      dplyr::n_distinct(dog_data_raw$file_key), "dogs.\n")
  
  # -- Step 4 : inner join with metadata --------------------------------------
  # Inner join: only dogs present in BOTH the CSV folder AND the metadata
  # are retained. Unmatched CSVs are reported as a warning.
  unmatched_csv <- setdiff(
    unique(dog_data_raw$join_key),
    unique(metadata$join_key)
  )
  
  if (length(unmatched_csv) > 0L) {
    warning(
      "load_dog_data(): ", length(unmatched_csv),
      " CSV file(s) have no metadata match and are EXCLUDED:\n",
      paste(sort(unmatched_csv), collapse = "\n"),
      "\nVerify filename conventions and metadata completeness."
    )
  }
  
  dog_data <- dog_data_raw |>
    dplyr::inner_join(metadata, by = "join_key")
  
  # Overwrite prefix_2 and setting with the values parsed from the filename
  # (the metadata join may introduce a duplicate prefix_2 from metadata cols)
  dog_data <- dog_data |>
    dplyr::mutate(
      prefix_2 = stringr::str_extract(file_key, "^[A-Z]{2}"),
      setting  = stringr::str_extract(file_key, "^[A-Z]{3}") |>
        stringr::str_sub(3, 3)
    )
  
  cat("load_dog_data(): after metadata join --",
      dplyr::n_distinct(dog_data$file_key), "dogs,",
      nrow(dog_data), "GPS fixes retained.\n")
  cat("load_dog_data(): CSV files with no metadata match (excluded):",
      length(unmatched_csv), "\n")
  
  dog_data
}