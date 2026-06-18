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
#' Function 2: Plot household sero-prevalence
#' 
#' @description This function uses individual GPS and sero-prevalence data to 
#' generate a sero-prevalence distribution layer through IDW interpolation
#' 
#' @param path Path used to save images
#' @param households_locations Data frame with three columns: id, long, lat.
#' @param status_df Data frame with columns: id, status (1 = positive, 0 = negative).
#' @param image_name Desired name of image for spatial seroprevalence distribution
#'
#' @export
household_seroprevalence <- function(path, 
                                     households_locations, 
                                     status_df, 
                                     image_name = "example_seroprevalence_distribution") {
  
  ### Initial variable/df checks
  
  # Check if columns of households_locations exist
  required_cols <- c("id", "long", "lat")
  missing_cols <- setdiff(required_cols, names(households_locations))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in households_locations:", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure coordinates are numeric
  if (!is.numeric(households_locations$lat) || !is.numeric(households_locations$lon)) 
    stop("Coordinates must be numeric.\n")
  
  # Check if status column in status_df exist
  if(!c("status") %in% names(status_df))
    stop("Missing sero-prevalence column 'status'.")
  
  # Check if status column is numeric
  if (!is.numeric(status_df$status)) 
    stop("Sero-prevalence must be numeric (binary).\n")

  
  ### Spatial distribution of household sero-prevalence
  
  # Merge movement with status
  full_df <- merge(household_locations, status_df, by = "id")
  
  # Turn into a spatial points dataframe
  sero_SPDF <- sp::SpatialPointsDataFrame(coords = full_df[,c("long", "lat")],
                                             data = full_df[,c("id", "status")])
  
  # Create object with only positives
  positives <- full_df[full_df$status==1,]
  
  # Create object with only negatives 
  negatives <- full_df[full_df$status==0,]
  
  # Set color
  color_sero <- leaflet::colorFactor(palette = c("blue", "red"), 
                                     domain = sero_SPDF$status)
  
  # Create leaflet map
  leaflet_map <- leaflet::leaflet(data = sero_SPDF) %>% 
    leaflet::addTiles() %>% 
    leaflet::addCircleMarkers(color = ~color_sero(status), 
                              stroke = FALSE,
                              fillOpacity = 0.8,
                              radius = 8)
  
  # Use mapview to create snapshot of leaflet map
  mapview::mapshot(
    leaflet_map, 
    file = file.path(path, paste0(image_name, ".png")),
    remove_controls = c("zoomControl", "layersControl"))
  
}
#'
#'
#'
#'
#' Function 3: Generate Seroprevalence Heatmap from Activity
#' 
#' @description
#' Function to generate seroprevalence heatmap from GPS and sero-data
#' 
#' @param movement_df Data frame of GPS trajectories (id, long, lat).
#' @param status_df Data frame of status (id, status), where status is 0/1.
#' @param base_raster A SpatRaster defining the study grid (e.g., from create_base_raster)
#' @param method either \code{relrisk} or inverse distance weighting (\code{IDW})
#' @param sigma Smoothing bandwidth (in degrees)
#' @param power either use power range (TRUE) or default power of 1 (FALSE)
#' 
#' @export
generate_seroprevalence_raster <- function(movement_df, 
                                           status_df, 
                                           base_raster,
                                           method, 
                                           sigma = 0.005,
                                           power = TRUE,
                                           path,
                                           image_name = "example_sero_map") {
  
  ### Initial variable/df checks
  
  # Check if columns of movement_df exist
  required_cols <- c("id", "long", "lat")
  missing_cols <- setdiff(required_cols, names(movement_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in movement_df:", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure coordinates are numeric
  if (!is.numeric(movement_df$lat) || !is.numeric(movement_df$lon)) 
    stop("Coordinates must be numeric.\n")
  
  # Check if status column in status_df exist
  if(!c("status") %in% names(status_df))
    stop("Missing sero-prevalence column 'status'.")
  
  # Check if status column is numeric
  if (!is.numeric(status_df$status)) 
    stop("Sero-prevalence must be numeric (binary).\n")
  
  # Check if base raster is raster layer
  if(class(base_raster) != "SpatRaster"){
    stop("'base_raster' has to be an object of class SpatRaster")
  }
  
  # Check the CRS and change to EPSG:4326 if necessary
  current_crs <- terra::crs(base_raster)
  if (current_crs == "") {
    message("Changing CRS to EPSG:4326")
    current_crs <- "EPSG:4326"
  }
  
  
  ### Prepare data for mapping
  
  # Join data
  full_df <- movement_df %>%
    dplyr::left_join(status_df, by = "id") %>%
    dplyr::mutate(status_fac = as.factor(as.character(status)))
  
  # Get SpatExtent of SpatRaster object
  ext_val <- as.vector(terra::ext(base_raster)) # xmin, xmax, ymin, ymax
  
  # Create observation window (owin)
  win <- spatstat.geom::owin(ext_val[1:2], ext_val[3:4])
  
  # Create marked point pattern object > cells with more GPS points have higher
  # weight in density calculation
  ppp_sero <- spatstat.geom::ppp(full_df$long, 
                               full_df$lat, 
                               window = win, 
                               marks = full_df$status)
  
  if(method == 'relrisk'){
    
    ### Risk mapping
    
    # Relative risk of being a case
    risk_map <- spatstat.explore::relrisk(X = ppp_sero, relative = TRUE)
    
    # Convert to spatraster
    sero_map <- terra::rast(risk_map)
    
    # Align with base raster
    sero_map <- terra::resample(sero_map, base_raster, method = "bilinear")
    
    # Safety clamp
    sero_map <- terra::clamp(sero_map, lower = 0, upper = 1)
    
    # Rename
    names(sero_map) <- "sero_probability"
    
  } 
  
  if(method == 'IDW'){
    
    
    ### Interpolation using Inverse Distance Weighting (IDW)
    
    # Overwrite power by calcuting power with lowest mean squared error
    if(power){
      
      # Set range of powers to test
      powers <- seq(0.05, 2, 0.05) 
      
      # Create empty vector of mean square error
      mse_vec <- NULL 
      
      # Compute IDW for multiple powers
      for(power in powers){idw_cv <- spatstat.explore::idw( X = ppp_sero, 
                                                            power=power, 
                                                            at="points",
                                                            se = FALSE) 
      
      # Save all mean squared errors
      mse_vec <- c(mse_vec, Metrics::mse(ppp_sero$marks,
                                idw_cv)) 
      }
      
      # Choose power with lowest mse
      idw_power <- powers[which.min(mse_vec)]
      
    } else {
      
      # If function argument power is FALSE -> set at 1
      idw_power <- 0.2
    }
    
    # Compute IDW
    sero_idw <- idw(ppp_sero, 
                    power = 0.2, 
                    at="pixels") 
    
    # Convert to spatraster
    sero_map <- terra::rast(sero_idw, crs = crs)
    
    # Align with base raster
    sero_map <- terra::resample(sero_map, base_raster, method = "bilinear")
    
    # Safety clamp
    sero_map <- terra::clamp(sero_map, lower = 0, upper = 1)
    
    # Rename
    names(sero_map) <- "sero_probability"
      
  }
  
  # Plot image
  png(paste0(file.path(path, image_name), ".png"))
  plot(sero_map)
  dev.off()
  
  # Return map
  return(sero_map)
}
#' 
#' 
#' Function 4: Create Integrated Environmental Exposure Heatmap
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
#' Function 5: Calculate Individual Dynamic Spatiotemporal Risk Index
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
  