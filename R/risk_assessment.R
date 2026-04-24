# =============================================================================
# R/environmental_risk.R
# Package : dogtrack
# Purpose : Quantify environmental risk using raster data
# Exported functions:
#   get_contact_probabilities()  -- uses GPS data to create contact probabilities
# Dependencies:
#   static.analysis.R
# =============================================================================
#' 
#' 
#' Function 1: Generate Seroprevalence Distribution Layer
#' 
#' @param movement_df Data frame with columns: id, long, lat.
#' @param status_df Data frame with columns: id, status (1 = positive, 0 = negative).
#' @param base_raster A SpatRaster defining the study grid.
#' @param sigma Smoothing parameter (Gaussian blur) for the intensity.
#'
#' @export
generate_seroprevalence_layer <- function(movement_df, status_df, base_raster, sigma = 1) {
  
  # Merge movement with status
  full_df <- merge(movement_df, status_df, by = "id")
  
  # Separate and rasterize
  pos_pts <- as.matrix(full_df[full_df$status == 1, c("long", "lat")])
  neg_pts <- as.matrix(full_df[full_df$status == 0, c("long", "lat")])
  
  pos_rast <- terra::rasterize(pos_pts, base_raster, fun = "count")
  neg_rast <- terra::rasterize(neg_pts, base_raster, fun = "count")
  
  # Fill zeros and smooth
  pos_rast <- terra::subst(pos_rast, NA, 0)
  neg_rast <- terra::subst(neg_rast, NA, 0)
  
  if(sigma > 0) {
    pos_rast <- terra::focal(pos_rast, w = 3, fun = "mean", na.rm = TRUE)
    neg_rast <- terra::focal(neg_rast, w = 3, fun = "mean", na.rm = TRUE)
  }
  
  # Calculate ratio (Positive distribution vs Total)
  sero_map <- pos_rast / (pos_rast + neg_rast + 1e-6)
  names(sero_map) <- "seroprevalence_score"
  return(sero_map)
}
#' 
#' 
#' 
#' Function 2: Create Integrated Environmental Exposure Heatmap
#' 
#' @param sero_layer Output from generate_seroprevalence_layer.
#' @param dog_density Raster of host density.
#' @param water_points sf object of water sources.
#' @param waste_points sf object of waste sources.
#' @param weights Named list: c(sero, density, env).
#'
#' @export
generate_environmental_exposure_map <- function(sero_layer, dog_density, 
                                                water_points = NULL, waste_points = NULL,
                                                weights = list(sero = 0.4, density = 0.3, env = 0.3)) {
  
  # Distance-based environmental risk
  env_risk <- terra::init(sero_layer, 0)
  if(!is.null(water_points)) env_risk <- env_risk + (1 / (1 + terra::distance(sero_layer, water_points)))
  if(!is.null(waste_points)) env_risk <- env_risk + (1 / (1 + terra::distance(sero_layer, waste_points)))
  
  # Helper to normalize 0-1
  norm <- function(x) (x - terra::minmax(x)[1]) / (terra::minmax(x)[2] - terra::minmax(x)[1])
  
  # Combined Score
  total_risk <- (weights$sero * norm(sero_layer)) + 
    (weights$density * norm(dog_density)) + 
    (weights$env * norm(env_risk))
  
  names(total_risk) <- "environmental_shedding_score"
  return(total_risk)
}
#'
#'
#' Function 3: Calculate Individual Dynamic Spatiotemporal Risk Index
#' 
#' @param total_locations Numeric. Number of unique cells visited in time window (e.g. 14 days).
#' @param shedding_score Numeric. Value from Function 2 at individual's current location.
#' @param weight_contact Numeric weight for social transmission.
#' @param weight_envir Numeric weight for environmental transmission.
#' @param loc_inf Binary (0/1). Does the current cell contain an infectious individual?
#' @param loc_sus Binary (0/1). Does the current cell contain a suspected individual?
#'
#' @export
dynamic_risk_index <- function(total_locations,
                               shedding_score,
                               weight_contact,
                               weight_envir,
                               loc_inf,
                               loc_sus) {
  
  # Wang et al. (2020) contact component
  # Ratio of risky locations to total locations visited
  contact_risk <- (loc_inf + loc_sus) / total_locations
  
  # Integrated formula
  # RI(x) = [Wc * (Contact Component)] + [We * (Environmental Component)]
  ri_x <- (weight_contact * contact_risk) + (weight_envir * shedding_score)
  
  return(ri_x)
}
#' create an infection hotspot layer from seroprevalence and GPS data
#'
#' @description Directly transforms a data frame of animal movements and status 
#' into a SpatRaster representing presence intensity.
#'
#' @param df A data.frame/data.table with columns: \code{long}, \code{lat}, \code{status}.
#' @param base_raster A \code{SpatRaster} to define the grid geometry.
#' @param type Character. Either \code{"ratio"} (positive/negative), 
#' \code{"density"} (total count), or \code{"both"}.
#' @param sigma Numeric. If > 0, applies a Gaussian blur (smoothing) to the raster.
#'
#' @return A \code{SpatRaster} with one or more layers.
#' @export
#' @import terra
generate_presence_raster <- function(df, base_raster, type = "ratio", sigma = 0) {
  
  # 1. Faster Coordinate Extraction
  # Extract coords as a matrix for terra
  coords <- as.matrix(df[, c("long", "lat")])
  
  # 2. Separate by status for the ratio
  pos_coords <- coords[df$status == 1, , drop = FALSE]
  neg_coords <- coords[df$status == 0, , drop = FALSE]
  
  # 3. Blazing fast rasterization
  # 'fun = sum' on a vector of 1s effectively counts points per cell
  pos_rast <- terra::rasterize(pos_coords, base_raster, fun = "count")
  neg_rast <- terra::rasterize(neg_coords, base_raster, fun = "count")
  
  # Fill NAs with 0
  pos_rast[is.na(pos_rast)] <- 0
  neg_rast[is.na(neg_rast)] <- 0
  
  # 4. Optional Smoothing (Gaussian Blur)
  # This turns discrete points into a continuous 'heatmap'
  if (sigma > 0) {
    # w must be odd
    fw <- max(3, floor(sigma * 3) %2% 2 + floor(sigma * 3)) 
    pos_rast <- terra::focal(pos_rast, w = fw, fun = "mean", na.rm = TRUE)
    neg_rast <- terra::focal(neg_rast, w = fw, fun = "mean", na.rm = TRUE)
  }
  
  # 5. Output Logic
  if (type == "ratio") {
    # Distribution of positive vs negative
    # Add small epsilon to avoid 0/0
    out <- pos_rast / (pos_rast + neg_rast + 1e-6)
    names(out) <- "pos_neg_distribution"
  } else if (type == "density") {
    out <- pos_rast + neg_rast
    names(out) <- "total_intensity"
  } else {
    out <- c(pos_rast, neg_rast)
    names(out) <- c("positive_intensity", "negative_intensity")
  }
  
  return(out)
}
#'
#' 
#' Function 2
#' Spatiotemporal Risk Heatmap that can be used in an individual-based disease 
#' transmission model. Here, we use different rasters with information about the 
#' environment and merge them to create a single layer that provides an estimate 
#' between 0 and 1 in the form of a heatmap on the degree of viral shedding as a 
#' proxy for environmental exposure to a pathogen. 
#'
#' @description Creates a normalized risk surface based on dog sero-prevalence, 
#' population density, and environmental point sources (water, waste, etc.).
#'
#' @param base_raster A \code{SpatRaster} template defining the extent and resolution.
#' @param dog_points An \code{sf} object containing dog locations. Must have a 
#' numeric \code{status} column (e.g., 0 for sero-negative, 1 for sero-positive).
#' @param dog_density Optional \code{SpatRaster}. If NULL, calculated via Kernel Density of \code{dog_points}.
#' @param water_points Optional \code{sf} points representing water sources.
#' @param waste_points Optional \code{sf} points representing waste/contamination sources.
#' @param weights A named list of weights: \code{dog}, \code{prevalence}, and \code{env}. Must sum to 1.
#'
#' @return A \code{SpatRaster} representing the integrated risk score (0 to 1).
#' @export
#' @import terra
#' @import sf
#' @importFrom gstat gstat
#' @importFrom stats predict
generate_risk_heatmap <- function(base_raster, 
                                  dog_points, 
                                  dog_density = NULL,
                                  water_points = NULL,
                                  waste_points = NULL,
                                  weights = list(dog = 0.4, env = 0.3, prevalence = 0.3)) {
  
  # 1. Validation -----------------------------------------------------------
  if (!inherits(base_raster, "SpatRaster")) stop("base_raster must be a SpatRaster.")
  if (!"status" %in% names(dog_points)) stop("dog_points must contain a 'status' column.")
  
  # Helper to normalize rasters to 0-1 scale
  .norm_rast <- function(x) {
    mm <- terra::minmax(x)
    if (mm[1] == mm[2]) return(x * 0)
    (x - mm[1]) / (mm[2] - mm[1])
  }
  
  # 2. Seroprevalence Interpolation -----------------------------------------
  # Turn point-based infection status into a continuous probability surface
  # Using Inverse Distance Weighting (IDW)
  gs <- gstat::gstat(formula = status ~ 1, locations = dog_points)
  prev_map <- terra::interpolate(base_raster, gs, debug.level = 0)
  prev_map <- .norm_rast(prev_map)
  
  # 3. Density Modeling -----------------------------------------------------
  if (is.null(dog_density)) {
    # Generate density if not provided (simple rasterization + Gaussian blur)
    # This is a 'package-safe' alternative to spatstat for density
    dog_rast <- terra::rasterize(dog_points, base_raster, fun = "count")
    dog_rast[is.na(dog_rast)] <- 0
    dog_density <- terra::focal(dog_rast, w = 5, fun = "mean", na.rm = TRUE)
  }
  dog_density <- .norm_rast(dog_density)
  
  # 4. Environmental Proximity ----------------------------------------------
  # We model risk as a 1/(1+d) decay from points
  env_risk <- terra::init(base_raster, 0)
  
  calc_dist_risk <- function(pts) {
    d <- terra::distance(base_raster, pts)
    return(1 / (1 + d))
  }
  
  if (!is.null(water_points)) {
    env_risk <- env_risk + calc_dist_risk(water_points)
  }
  if (!is.null(waste_points)) {
    env_risk <- env_risk + calc_dist_risk(waste_points)
  }
  
  env_risk <- .norm_rast(env_risk)
  
  # 5. Final Calculation ----------------------------------------------------
  risk_score <- (weights$prevalence * prev_map) + 
    (weights$dog * dog_density) + 
    (weights$env * env_risk)
  
  names(risk_score) <- "integrated_risk_index"
  return(risk_score)
}
#'
#'
#' Function 3
#' Calculate exposure risk for each individual based on spatiotemporal GPS data
#' as well as interaction information. The risk index is based on Wang et al. (2020)
#' (http://mhealth.jmir.org/2020/6/e19457/) and is quantified as the ratio of the
#' weighted proportion of locations containing infectious and suspected infectious
#' individuals. 
#' 
#' RI(x) = the individual dynamic spatiotemporal risk index for individual x. Lxi 
#' and Lxs are the locations that individual x has visited in the past 14 days 
#' which contain infectious or suspected individuals. Lxi:Lxs / sum(L)t then gives 
#' the ratio of these locations to all locations that have been visited by 
#' individual x. We complement this RI(x) with 
#' This function uses latitude and longitude coordinates and an environmental risk
#' map as well as contact probabilities to establish a risk probility for contagion.
#'
#' @param total_locations the total number of unique raster cells visited by x
#' @param shedding_score Represents viral load in environment of location 1 (0-1)
#' @param weight_contact a weighting factor for direct contact risk
#' @param weight_envir a weighting factor for environmental risk
#' @param loc_inf location contains infectious individuals (0/1)
#' @param loc_sus location contains suspected infectious individuals (0/1)
#' 
#' @include contact_probabilities.R
#'
#' @return A tibble with columns `dog_a`, `dog_b`, `t_a`, `t_b`, `dist`, `p_contact`.
#' 
#' @export
dynamic_risk_index <- function(total_locations,
                               shedding_score,
                               weight_contact,
                               weight_envir,
                               loc_inf,
                               loc_sus) {
  