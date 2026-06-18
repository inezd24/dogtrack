# =============================================================================
# R/metadata.R
# Package : dogtrack
# Purpose : Read metadata from a workbook or flat file and build
#           the prefix-to-site routing table used by the pipeline.
#
# Exported functions:
#   load_metadata()    -- reads metadata, normalises columns, builds join_key
#   build_prefix_map() -- derives prefix_2 -> field_site lookup
#
# Internal helpers (not exported):
#   .find_header_row()        -- XLSX/XLS only: locates header row above title rows
#   .harmonise_column_names() -- maps built-in language/abbreviation variants
#   .harmonise_site_values()  -- normalises field_site spelling variants
#   .read_and_harmonise()     -- snake_cases names, applies variants + col_aliases
#   .clean_rows()             -- drops incomplete rows, builds join_key/prefix_2
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported -- XLSX/XLS only
# Scans the first `max_rows` rows of an Excel sheet to find the row containing
# the actual column headers. Detection is based on a regex that matches common
# ways of labelling a household/housing unit across languages (English, French,
# Spanish, Swahili, etc.). Returns the 1-based row number for read_excel(skip=).
# The `target` regex can be overridden via `header_target` in load_metadata()
# if none of the built-in patterns match your workbook.
# Not called for .csv / .txt files, which have no title rows above the header.
# -----------------------------------------------------------------------------
.find_header_row <- function(metadata_path, sheet,
                             target = paste(
                               "household",  # English
                               "maison",     # French
                               "hogar",      # Spanish
                               "kaya",       # Swahili
                               "hh[^a-z]",   # common abbreviation: HH_ID, HH No
                               "house",      # informal English
                               sep = "|"
                             ),
                             max_rows = 20L) {
  
  raw <- readxl::read_excel(
    metadata_path, sheet     = sheet,
    col_names = FALSE, n_max = max_rows,
    col_types = "text"
  )
  
  hit <- which(apply(raw, 1, function(row) {
    any(stringr::str_detect(
      tidyr::replace_na(row, ""),
      stringr::regex(target, ignore_case = TRUE)
    ))
  }))
  
  if (length(hit) == 0L) {
    stop(
      "load_metadata(): cannot locate a header row in sheet '", sheet,
      "' (searched first ", max_rows, " rows).\n",
      "The default search looks for terms like 'household', 'maison', 'hogar', ",
      "'kaya', 'HH', or 'house'.\n",
      "If your workbook uses a different label, pass `header_target` to ",
      "load_metadata(), e.g.: header_target = 'respondent_id'"
    )
  }
  
  hit[[1L]]
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Harmonises known column NAME variants (e.g. French/Swahili abbreviations)
# to the standard canonical names used by the rest of the pipeline. Only
# renames a variant if the canonical name is NOT already present. This runs
# BEFORE col_aliases and BEFORE the required-column check, so that built-in
# variants are resolved automatically without requiring the user to specify
# col_aliases for languages/conventions already known to the package.
# To support a new language or naming convention, add entries below -- no
# other changes are needed.
# -----------------------------------------------------------------------------
.harmonise_column_names <- function(df) {
  
  variants <- list(
    household_id = c("maison_id", "hh_id", "hh_no", "house_id"),
    field_site   = c("ville", "site", "location", "study_site"),
    dog_id       = c("nr_chien", "no_chien", "dog_number", "dog_no")
  )
  
  for (standard in names(variants)) {
    if (!standard %in% names(df)) {
      found <- intersect(variants[[standard]], names(df))
      if (length(found) > 0L) {
        df <- dplyr::rename(df, !!standard := !!found[[1L]])
      }
    }
  }
  
  df
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Normalises known field_site spelling/casing variants to their canonical
# form (e.g. "ndjamena" -> "N'Djamena"). Operates on VALUES, not column names;
# assumes `field_site` already exists. Add new entries here whenever a new
# site spelling variant is encountered.
# -----------------------------------------------------------------------------
.harmonise_site_values <- function(df) {
  
  site_canonical <- c(
    "n'djamena" = "N'Djamena",
    "n djamena" = "N'Djamena",
    "ndjamena"  = "N'Djamena",
    "masaka"    = "Masaka",
    "arua"      = "Arua",
    "soroti"    = "Soroti"
  )
  
  if ("field_site" %in% names(df)) {
    df <- df |>
      dplyr::mutate(
        field_site = {
          key   <- stringr::str_to_lower(stringr::str_trim(field_site))
          canon <- site_canonical[key]
          dplyr::if_else(!is.na(canon), canon, field_site)
        }
      )
  }
  
  df
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Resolves a raw data frame's column names to the canonical names expected
# by the pipeline, in this order:
#   1. snake_case normalisation (lower, non-alphanumeric -> "_", collapse,
#      strip leading/trailing underscores)
#   2. built-in language/abbreviation variants, via .harmonise_column_names()
#   3. caller-supplied `col_aliases` (a named list mapping canonical names to
#      whatever the file actually uses) -- aliases are themselves
#      snake_case-normalised before matching
# The three canonical names this function resolves toward are:
#   "household_id", "dog_id", "field_site"
# A column already named canonically at any step is left untouched.
# `col_aliases` is expected to be validated by the caller (load_metadata).
# -----------------------------------------------------------------------------
.read_and_harmonise <- function(df_raw, source_label, col_aliases,
                                required_canonical) {
  
  .to_snake <- function(x) {
    stringr::str_remove(
      stringr::str_replace_all(
        stringr::str_replace_all(stringr::str_to_lower(x), "[^a-z0-9]", "_"),
        "_+", "_"),
      "^_|_$")
  }
  
  # Step 1: snake_case
  df <- dplyr::rename_with(
    dplyr::select(df_raw, dplyr::where(~!all(is.na(.)))),
    .to_snake
  )
  
  # Step 2: built-in language/abbreviation variants
  df <- .harmonise_column_names(df)
  
  # Step 3: user-supplied col_aliases
  if (!is.null(col_aliases)) {
    for (canonical in names(col_aliases)) {
      alias_snake <- .to_snake(col_aliases[[canonical]])
      if (alias_snake %in% names(df) && !canonical %in% names(df)) {
        df <- dplyr::rename(df, !!canonical := !!alias_snake)
      }
    }
  }
  
  df
}


# -----------------------------------------------------------------------------
# Internal helper -- NOT exported
# Drops rows missing any required column, trims whitespace, and builds the
# derived `join_key` and `prefix_2` columns. By this point column NAMES are
# already fully resolved (snake_case + built-in variants + col_aliases), so
# only value-level site-name canonicalisation remains.
# -----------------------------------------------------------------------------
.clean_rows <- function(df, source_label) {
  df <- .harmonise_site_values(df)   # value-level only; names already resolved
  n_raw <- nrow(df)
  df <- dplyr::mutate(
    dplyr::filter(df, !is.na(household_id), !is.na(dog_id), !is.na(field_site)),
    household_id = stringr::str_trim(as.character(household_id)),
    dog_id       = stringr::str_trim(as.character(dog_id)),
    field_site   = stringr::str_trim(as.character(field_site)),
    sheet_name   = source_label,
    join_key     = dplyr::if_else(
      stringr::str_detect(dog_id, "-"),
      dog_id,
      paste0(household_id, "-",
             stringr::str_pad(stringr::str_remove(dog_id, "\\.0$"),
                              width = 2, side = "left", pad = "0"))
    ),
    prefix_2 = stringr::str_extract(household_id, "^[A-Z]{2}")
  )
  n_dropped <- n_raw - nrow(df)
  if (n_dropped > 0) {
    warning("load_metadata(): '", source_label, "': dropped ", n_dropped,
            " row(s) with missing household_id, dog_id, or field_site.")
  }
  df
}


# -----------------------------------------------------------------------------
#' Load and harmonise metadata from a workbook or flat file
#'
#' Reads metadata from an Excel workbook (`.xlsx`, `.xls`), CSV, or
#' tab-delimited text file. For Excel workbooks, all sheets containing the
#' three required columns are loaded and row-bound. Column names are normalised
#' to `snake_case` immediately after reading; known language variants (French,
#' Swahili, etc.) are mapped to canonical names automatically.
#'
#' @param metadata_path Character. Path to the metadata file.
#'   Supported formats: `.csv`, `.txt` (tab-delimited), `.xlsx`, `.xls`.
#' @param col_aliases Named list, or `NULL`. Use this when your file's column
#'   names differ from the canonical names expected by the pipeline
#'   (`household_id`, `dog_id`, `field_site`) AND are not already covered by
#'   the package's built-in variants (see `.harmonise_column_names()`). Map
#'   each canonical name to the name used in your file, e.g.:
#'   ```r
#'   col_aliases = list(
#'     household_id = "HH_ID",
#'     dog_id       = "DogNo",
#'     field_site   = "Site"
#'   )
#'   ```
#'   You only need to supply entries for columns that differ; omit any that
#'   already match or are already covered by a built-in variant. Aliases are
#'   matched after snake_case normalisation, so punctuation and capitalisation
#'   in the alias value do not matter (`"Dog No."`, `"dog_no"`, and `"DOG NO"`
#'   all resolve identically).
#' @param sheet Character or integer, or `NULL`. Excel files only. Pin a
#'   specific sheet by name or position, bypassing auto-detection. Useful once
#'   you know your file structure and want to avoid the scanning overhead.
#' @param header_target Character. Excel files only. Regex used to locate the
#'   header row in each sheet (see `.find_header_row()`). The default covers
#'   common words for household across several languages (English, French,
#'   Spanish, Swahili) and the abbreviation `HH`. Override if your workbook
#'   uses none of these, e.g. `header_target = "respondent_id"`.
#'
#' @return A tibble with one row per dog deployment. The following columns are
#'   always present in the output, using the canonical names listed below
#'   regardless of how they are named in the source file. If your file uses
#'   names not already covered by the package's built-in variants, you must
#'   supply them via `col_aliases` -- the function will error if any required
#'   column cannot be resolved:
#'   \describe{
#'     \item{`household_id`}{Character. Unique household identifier.}
#'     \item{`dog_id`}{Character. Dog identifier within the household.}
#'     \item{`field_site`}{Character. Canonical site name (e.g. `"Masaka"`,
#'       `"N'Djamena"`). Known spelling variants are normalised automatically;
#'       unrecognised values are kept as-is.}
#'     \item{`join_key`}{Character. Constructed as
#'       `"{household_id}-{zero-padded dog_id}"`. Used to match GPS CSV
#'       filenames to metadata rows. If `dog_id` already contains a hyphen it
#'       is used directly.}
#'     \item{`prefix_2`}{Character. First two uppercase letters of
#'       `household_id` (e.g. `"UM"` for Uganda Masaka). Used for file routing
#'       by [build_prefix_map()].}
#'     \item{`sheet_name`}{Character. Source sheet (Excel) or filename
#'       (CSV/TXT), retained for traceability.}
#'   }
#'   Additional columns present in the source file are carried through
#'   unchanged.
#'
#' @details
#' **Column name resolution** happens in three steps, in order:
#' 1. All column names are converted to `snake_case`.
#' 2. Known language/abbreviation variants (French, Swahili, common
#'    abbreviations like `HH_ID`) defined in `.harmonise_column_names()` are
#'    mapped to canonical names automatically.
#' 3. Any remaining non-standard names must be declared via `col_aliases`.
#'
#' If after all three steps a required column is still absent, the function
#' stops with an informative error listing what was found and how to fix it.
#'
#' **Sheet selection (Excel only):** every sheet is scanned; those containing
#' all three required columns (after name resolution) are loaded and row-bound.
#' Sheets lacking any required column are silently skipped. Use `sheet` to pin
#' a specific sheet and skip scanning entirely.
#'
#' **Header detection (Excel only):** workbooks often have one or more title
#' rows above the actual column headers. `.find_header_row()` locates the
#' header row by matching a regex against cell values; the default pattern
#' covers common household-related terms across languages. Override with
#' `header_target` if needed. This step does not apply to CSV/TXT files,
#' where the first row is always treated as the header.
#'
#' Rows missing any of `household_id`, `dog_id`, or `field_site` after
#' cleaning are dropped with a warning.
#'
#' @examples
#' \dontrun{
#' # Standard Excel workbook -- auto-detects qualifying sheets
#' meta <- load_metadata("data/metadata.xlsx")
#'
#' # File with non-standard column names
#' meta <- load_metadata(
#'   "data/metadata.xlsx",
#'   col_aliases = list(household_id = "HH_ID", dog_id = "DogNo")
#' )
#'
#' # CSV input
#' meta <- load_metadata("data/metadata.csv")
#'
#' dplyr::count(meta, field_site)
#' }
#'
#' @export
load_metadata <- function(
    metadata_path,
    col_aliases   = NULL,
    sheet         = NULL,
    header_target = paste(
      "household", "maison", "hogar", "kaya", "hh[^a-z]", "house",
      sep = "|"
    )
) {
  
  # ── 0. File existence and type check ────────────────────────────────────────
  if (!fs::file_exists(metadata_path)) {
    stop("load_metadata(): file not found: ", metadata_path)
  }
  
  ext <- tolower(tools::file_ext(metadata_path))
  if (!ext %in% c("csv", "txt", "xlsx", "xls")) {
    stop(
      "load_metadata(): unsupported file type '.", ext, "'.\n",
      "Supported types: .csv, .txt, .xlsx, .xls"
    )
  }
  
  required_canonical <- c("household_id", "dog_id", "field_site")
  
  # ── 1. Validate col_aliases once, up front ──────────────────────────────────
  if (!is.null(col_aliases)) {
    if (!is.list(col_aliases) || is.null(names(col_aliases))) {
      stop(
        "load_metadata(): `col_aliases` must be a named list, e.g.\n",
        "  list(household_id = 'HH_ID', dog_id = 'DogNo')"
      )
    }
    bad_keys <- setdiff(names(col_aliases), required_canonical)
    if (length(bad_keys) > 0) {
      warning(
        "load_metadata(): `col_aliases` contains unrecognised key(s): ",
        paste(bad_keys, collapse = ", "),
        ".\nExpected keys: ", paste(required_canonical, collapse = ", ")
      )
    }
  }
  
  # ── 2. Helper: validate required columns present ────────────────────────────
  .check_required <- function(df, source_label) {
    missing <- setdiff(required_canonical, names(df))
    if (length(missing) > 0) {
      stop(
        "load_metadata(): '", source_label, "' is missing column(s): ",
        paste(missing, collapse = ", "), ".\n",
        "Columns found (after harmonisation): ", paste(names(df), collapse = ", "), ".\n",
        "If your columns have non-standard names, supply `col_aliases`, e.g.:\n",
        "  col_aliases = list(dog_id = 'DogNo', field_site = 'Site')"
      )
    }
    invisible(df)
  }
  
  # ── 3. Branch on file type ───────────────────────────────────────────────────
  
  # ── 3a. CSV / TXT ────────────────────────────────────────────────────────────
  if (ext %in% c("csv", "txt")) {
    sep   <- if (ext == "csv") "," else "\t"
    label <- fs::path_file(metadata_path)
    df_raw <- suppressMessages(
      readr::read_delim(metadata_path, delim = sep,
                        col_types = readr::cols(.default = "c"),
                        show_col_types = FALSE)
    )
    df <- .read_and_harmonise(df_raw, source_label = label,
                              col_aliases = col_aliases,
                              required_canonical = required_canonical)
    .check_required(df, source_label = label)
    cat("load_metadata(): reading flat file '", label, "'\n", sep = "")
    return(.clean_rows(df, source_label = label))
  }
  
  # ── 3b. XLSX / XLS ───────────────────────────────────────────────────────────
  sheets_all <- readxl::excel_sheets(metadata_path)
  if (length(sheets_all) == 0) {
    stop("load_metadata(): the workbook at '", metadata_path, "' has no sheets.")
  }
  
  if (!is.null(sheet)) {
    # Caller pinned a specific sheet -- use it directly, no scanning
    target_sheets <- if (is.numeric(sheet)) sheets_all[[sheet]] else sheet
    if (!target_sheets %in% sheets_all) {
      stop("load_metadata(): sheet '", target_sheets, "' not found in '",
           metadata_path, "'.\n",
           "Sheets present: ", paste(sheets_all, collapse = ", "))
    }
    cat("load_metadata(): using pinned sheet '", target_sheets, "'\n", sep = "")
  } else {
    # Auto-detect: scan headers of every sheet, keep those with all three columns
    target_sheets <- purrr::keep(sheets_all, function(s) {
      tryCatch({
        header_row <- .find_header_row(metadata_path, s, target = header_target)
        df_raw <- suppressMessages(
          readxl::read_excel(metadata_path, sheet = s,
                             skip = header_row - 1L, col_types = "text", n_max = 0)
        )
        df_peek <- .read_and_harmonise(df_raw, source_label = s,
                                       col_aliases = col_aliases,
                                       required_canonical = required_canonical)
        all(required_canonical %in% names(df_peek))
      }, error = function(e) FALSE)
    })
    
    if (length(target_sheets) == 0L) {
      stop(
        "load_metadata(): no sheet in '", metadata_path,
        "' contains all three required columns\n",
        "  (", paste(required_canonical, collapse = ", "), ") ",
        "after name harmonisation.\n",
        "Sheets scanned: ", paste(sheets_all, collapse = ", "), ".\n",
        "If your columns have non-standard names, supply `col_aliases`, e.g.:\n",
        "  col_aliases = list(household_id = 'HH_ID', dog_id = 'DogNo')\n",
        "If your header row uses a different label, supply `header_target`, e.g.:\n",
        "  header_target = 'respondent_id'"
      )
    }
    
    cat("load_metadata(): auto-detected", length(target_sheets),
        "qualifying sheet(s):", paste(target_sheets, collapse = ", "), "\n")
  }
  
  # Full read of all qualifying sheets, row-bound into one tibble
  purrr::map_dfr(target_sheets, function(s) {
    header_row <- .find_header_row(metadata_path, s, target = header_target)
    df_raw <- suppressMessages(
      readxl::read_excel(metadata_path, sheet = s,
                         skip = header_row - 1L, col_types = "text")
    )
    df <- .read_and_harmonise(df_raw, source_label = s,
                              col_aliases = col_aliases,
                              required_canonical = required_canonical)
    .check_required(df, source_label = s)
    .clean_rows(df, source_label = s)
  })
  
}  # end load_metadata()


# -----------------------------------------------------------------------------
#' Build a prefix-to-site routing table from loaded metadata
#'
#' Derives a lookup table mapping each 2-character filename prefix (e.g. `"UM"`,
#' `"TN"`) to its `field_site` (e.g. `"Masaka"`, `"N'Djamena"`). This table
#' is used internally by the pipeline to route GPS CSV files to their correct
#' site-specific thresholds without any hardcoding.
#'
#' @param metadata A tibble as returned by [load_metadata()], containing at
#'   minimum the columns `prefix_2` and `field_site`.
#'
#' @return A tibble with columns `prefix_2` and `field_site`. One row per
#'   unique prefix.
#'
#' @details
#' Stops with an informative error if any 2-character prefix maps to more than
#' one `field_site` -- this would make routing ambiguous and indicates that two
#' sites share a country + city code, which should not happen.
#'
#' @examples
#' \dontrun{
#' meta    <- load_metadata("data/metadata.xlsx")
#' pfx_map <- build_prefix_map(meta)
#' # prefix_2   field_site
#' # "UM"        "Masaka"
#' # "UA"        "Arua"
#' # "TN"        "N'Djamena"
#' }
#'
#' @export
build_prefix_map <- function(metadata) {
  
  map <- metadata |>
    dplyr::distinct(prefix_2, field_site) |>
    dplyr::filter(!is.na(prefix_2))
  
  ambiguous <- map |>
    dplyr::count(prefix_2) |>
    dplyr::filter(n > 1L) |>
    dplyr::pull(prefix_2)
  
  if (length(ambiguous) > 0) {
    stop(
      "build_prefix_map(): the following prefix(es) map to more than one ",
      "field_site -- routing is ambiguous:\n",
      paste(ambiguous, collapse = ", "), "\n",
      "Ensure that each site has a unique country + city code combination."
    )
  }
  
  map
}