# =============================================================================
# R/calculate_contact_prob.R
# Package : dogtrack
# Purpose : High-performance calculation of contact probabilities using 
#           spatiotemporal decay and optional path-awareness.
# Exported functions:
#   get_contact_probabilities()  -- uses GPS data to create contact probabilities
# Dependencies:
#   static.analysis.R
# =============================================================================

#' Calculate contact probabilities with aggregation and path awareness
#'
#' This function identifies potential contacts between individuals using a
#' forward-looking temporal window. It supports high-frequency data aggregation,
#' and path-aware distance checks.
#'
#' @param dog_data A tibble or data frame containing dog GPS fixes.
#' @param id_col Character. Name of the ID column. Default `"file_key"`.
#' @param time_col Character. Name of the datetime column. Default `"datetime"`.
#' @param lat_col Numeric. Name of the latitude/Y column. Default `"lat"`.
#' @param lon_col Numeric. Name of the longitude/X column. Default `"lon"`.
#' @param field_site_col Character. Name of site column. Required if `site_specific = TRUE`.
#' @param site_specific Logical. If `TRUE`, matches `r` to the `field_site`.
#' @param r Either a numeric global buffer or a tibble from [derive_site_thresholds()].
#' @param time_lag Logical. If `TRUE`, allows asynchronous contacts within `delta_t_max`.
#' @param agg_window Numeric. Seconds to aggregate pings. Only used if `time_lag = FALSE`.
#' @param delta_t_max Numeric. Maximum time lag (seconds) for contact. Default `30`.
#' @param v_max Numeric. Max assumed velocity (m/s). Default `5.6` (approx 20km/h).
#' @param lambda Numeric. Decay constant. If `NULL` and `path_aware = TRUE`, becomes adaptive.
#' @param path_aware Logical. If `TRUE`, uses a 3-point window (prev, curr, next) for distances.
#' 
#' @include static_analysis.R
#'
#' @return A tibble with columns `dog_a`, `dog_b`, `t_a`, `t_b`, `dist`, `p_contact`.
#' 
#' @export
get_contact_probabilities <- function(dog_data, 
                                      id_col = "file_key", 
                                      time_col = "datetime", 
                                      lat_col = "lat", 
                                      lon_col = "lon",
                                      field_site_col = "field_site",
                                      site_specific = FALSE,
                                      r = 5, 
                                      time_lag = TRUE,
                                      agg_window = 10,
                                      delta_t_max = 30, 
                                      v_max = 5.6, 
                                      lambda = 2.0,
                                      path_aware = FALSE) {
  
  # 1. Internal Validation and Column Normalization ----------------------------
  # ----------------------------------------------------------------------------
  
  # Check if columns exist
  required_cols <- c(id_col, time_col, lat_col, lon_col)
  if (site_specific) required_cols <- c(required_cols, field_site_col)
  
  # Document missing columns
  missing_cols <- setdiff(required_cols, names(dog_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in dog_data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Subset and Rename
  dt <- data.table::as.data.table(dog_data)[, ..required_cols]
  
  # Rename for consistency
  if (site_specific) {
    data.table::setnames(dt, old = required_cols, new = c("id", "datetime", "lat", "lon", "field_site"))
  } else {
    data.table::setnames(dt, old = required_cols, new = c("id", "datetime", "lat", "lon"))
  }
  
  # Format check
  if (!is.numeric(dt$lat) || !is.numeric(dt$lon)) 
    stop("Coordinates must be numeric.\n")
  
  # Format check
  if (!lubridate::is.POSIXct(dt$datetime)) {
    stop("Time column must be a POSIXct object (e.g., '2023-10-25 14:30:05').\n")
  }
  
  # Ensure data isn't just dates (must have H:M:S)
  if (all(format(dt$datetime, "%H%M%S") == "000000")) {
    warning("Timestamps appear to be dates only (midnight). H:M:S is required for contact calculation.\n")
  }
  
  # Check and report defaults
  if ((time_lag) && is.null(agg_window)){
    message("Time lag is TRUE but no aggregate window provided. Using default of 10 seconds.\n")
  }
  if (is.null(delta_t_max)){
    message("Maximum time window (temporal) lag is not provided. Default of 30 seconds used.\n")
  }
  if (is.null(r)){
    message("Spatial error is not provided. Default of 5 meters is used.\n")
  }
  r_label <- if(site_specific) "Empirical (Site-specific)" else paste0(r, "m (Global)")
  if (is.null(v_max)){
    message("Maximum velocity is not provided. Default of 5.6 meters per second is used.\n")
  }
  if (is.null(lambda) && !path_aware) {
    message("No lambda provided and path_aware is FALSE. Using default lambda = 2.0.\n")
  }
  if (is.null(lambda)) {
    lambda_label <- if(path_aware) "Adaptive (Movement-based)" else {
      "2.0 (Default fallback)"
    }
  } else {
    lambda_label <- as.character(lambda)
  }
  
  # Give summary of values at start of function:
  param_msg <- paste0(
    "Running contact analysis:\n",
    "  - Buffer (r): ", r_label, "\n",
    "  - Max Lag (dt): ", delta_t_max, "s\n",
    "  - Max Velocity (v): ", v_max, "m/s\n",
    "  - Lambda: ", lambda_label, "\n",
    "  - Path Aware: ", path_aware, "\n",
    "  - Time Lag Mode: ", time_lag, if(!time_lag) paste0(" (Agg: ", agg_window, "s)") else "", "\n"
  )
  message(param_msg)
  
  
  # 2. Pre-processing: Aggregation ---------------------------------------------
  # ----------------------------------------------------------------------------
  
  # Only if there is no time-lag and aggregation window was set
  if (!time_lag && agg_window > 0) {
    # 1. Round to absolute time buckets (e.g., every 10s from 00:00:00)
    # This ensures Dog A at 09:00:02 and Dog B at 09:00:08 both become 09:00:10
    dt[, datetime := as.POSIXct(ceiling(as.numeric(datetime) / agg_window) * agg_window, 
                                origin = "1970-01-01")]
    
    # 2. Extract Date to ensure we don't aggregate across midnight
    dt[, date_only := as.Date(datetime)]
    
    # 3. Aggregate by Individual AND Date AND Time-bucket
    dt <- dt[, .(lat = mean(lat, na.rm = TRUE), 
                 lon = mean(lon, na.rm = TRUE)), 
             by = .(id, date_only, datetime)]
    
    # Cleanup temp column
    dt[, date_only := NULL]
  }
  
  # 3. Path awareness ----------------------------------------------------------
  # ----------------------------------------------------------------------------
  if (path_aware) {
    data.table::setkey(dt, id, datetime)
    dt[, `:=`(
      prev_lat = data.table::shift(lat, type = "lag"), 
      prev_lon = data.table::shift(lon, type = "lag"),
      next_lat = data.table::shift(lat, type = "lead"),
      next_lon = data.table::shift(lon, type = "lead")
    ), by = id]
  }
  
  # 4. Spatiotemporal Join -----------------------------------------------------
  # ----------------------------------------------------------------------------
  # Prepare time keys for non-equi join
  dt[, t_search := datetime]
  
  if (time_lag) {
    # Look forward from 'datetime' to 'datetime + delta_t_max'
    dt[, t_end := datetime + delta_t_max]
    contacts <- dt[dt, on = .(t_search >= datetime, t_search <= t_end), 
                   allow.cartesian = TRUE, nomatch = 0]
  } else {
    # Synchronized exact match
    contacts <- dt[dt, on = .(datetime), allow.cartesian = TRUE, nomatch = 0]
  }
  
  # Remove self-contact and duplicate dyads (A-B and B-A)
  contacts <- contacts[id < i.id]
  
  # 5. Distance Calculation & Path Logic ---------------------------------------
  # ----------------------------------------------------------------------------
  # Calculate Current-to-Current (Haversine) distance
  contacts[, dist_curr := sqrt((lat - i.lat)^2 + (lon - i.lon)^2)]
  
  if (path_aware) {
    # Calculate the cross-distances for the 3-point window
    contacts[, dist_prev := sqrt((prev_lat - i.prev_lat)^2 + (prev_lon - i.prev_lon)^2)]
    contacts[, dist_next := sqrt((next_lat - i.next_lat)^2 + (next_lon - i.next_lon)^2)]
    
    # Determine the minimum distance encountered across the "step"
    contacts[, dist_min := pmin(dist_curr, dist_prev, dist_next, 
                                sqrt((prev_lat - i.lat)^2 + (prev_lon - i.lon)^2),
                                sqrt((lat - i.prev_lat)^2 + (lon - i.prev_lon)^2),
                                na.rm = TRUE)]
    
    # ADAPTIVE LAMBDA CALCULATION
    # If the user didn't provide a manual lambda, we calculate it based on trend
    if (is.null(lambda)) {
      contacts[, lambda_val := data.table::fcase(
        dist_next < dist_curr, 1.5,   # Converging: More forgiving (Lower lambda)
        dist_next > dist_curr, 3.0,   # Diverging: Stricter (Higher lambda)
        default = 2.0                 # Neutral/Static
      )]
    } else {
      contacts[, lambda_val := lambda] # Use manual user input
    }
  } else {
    contacts[, dist_min := dist_curr]
    contacts[, lambda_val := if(is.null(lambda)) 2.0 else lambda]
  }
  
  # 6. Setting spatial error ---------------------------------------------------
  # ----------------------------------------------------------------------------
  if (site_specific) {
    if (is.numeric(r)) stop("site_specific is TRUE but r is numeric. Provide thresholds tibble.")
    thresh_dt <- data.table::as.data.table(r)[, .(field_site, r_val = dist_noise_m)]
    dt <- merge(dt, thresh_dt, by = "field_site", all.x = TRUE)
    dt[is.na(r_val), r_val := 5] # Fallback
  } else {
    dt[, r_val := as.numeric(r)]
  }
  
  # 7. Apply Probability Formula -----------------------------------------------
  # ----------------------------------------------------------------------------
  two_r <- 2 * r
  max_reach <- v_max * delta_t_max
  
  # Avoid division by zero if delta_t_max is 0
  if (max_reach == 0) max_reach <- 1e-6 
  
  contacts[, p_contact := data.table::fifelse(
    dist_min <= two_r, 1.0, 
    data.table::fifelse(
      dist_min > (two_r + max_reach), 0.0,
      exp(-lambda_val * (dist_min - two_r) / max_reach)
    )
  )]
  
  # 8. Final Cleanup and Return ------------------------------------------------
  # ----------------------------------------------------------------------------
  res <- contacts[, .(
    dog_a = id, 
    dog_b = i.id, 
    t_a = datetime, 
    t_b = i.datetime, 
    dist = dist_final, 
    p_contact = p_contact
  )]
  
  return(dplyr::as_tibble(res))
}