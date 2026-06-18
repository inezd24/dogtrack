#' Turn spatial data into rasters with integrated features
#' 
#' @description This function creates a standardized Base Raster from any spatial 
#' data with the option of including a single feature column and turns them into 
#' a raster using terra. You can provide the resolution, a buffer for the bounding 
#' box, the crs, and also decide how empty cells should be treated. Moreover, if 
#' there are multiple GLS points per cell, you can choose how terra should 
#' transform the corresponding values to a single value. If data is categorical, 
#' this is not possible.
#' 
#' @param spatial_data A data frame or an \code{sf} object. If data frame, it
#' should contain the following columns:
#' - 'long': The longitude of the GPS point (georeferenced)
#' - 'lat': The latitude of the GPS point (georeferenced)
#' The data frame can optionally also have a feature column with any column name. 
#' @param res Numeric. Resolution in decimal degrees (e.g., 0.001 ~ 110m).
#' @param buffer Numeric. Buffer around the bounding box (in degrees) to avoid 
#' edge effects.
#' @param crs Character. The CRS string, default is EPSG:4326.
#' @param feature_col Character. Name of an optional column with numeric or 
#' categorical values representing cell features
#' @param fun Function to use, e.g. mean, when multiple numeric values per cell
#' @param na_to_zero logical. Turn empty cells into zeros or leave as NA. 
#' 
#' @return A SpatRaster object
#'
#' @export
#' 
#' @examplesIf rlang::is_interactive()
#' # Get example spatial data from roads in Ndjamena (not corresponding to reality)
#' df_ndjamena <- read.csv(system.file("extdata/ndjamena_roads.csv", 
#'                                     package = "dogtrack"
#'                                     ), 
#'                         header = TRUE, 
#'                         row.names = 1, 
#'                         stringsAsFactors = FALSE
#')
#' # Turn spatial data into a raster with points
#' base_raster <- create_base_raster(
#'   spatial_data = df_ndjamena,
#'   res = 0.001, # default resolution 
#'   buffer = 0.02, # Add slightly bigger buffer than normal
#'   crs = "EPSG:4326",
#'   fun = 'max', # take the maximum value of all cell values
#'   na_to_zero = TRUE # turn NAs into zeros
#' )
create_base_raster <- function(spatial_data, 
                               res = 0.001, 
                               buffer = 0.01, 
                               crs = "EPSG:4326",
                               feature_col = NULL,
                               fun = 'mean',
                               na_to_zero = FALSE) {
  
  # If data is not sf object, turn into sf object
  if (!inherits(spatial_data, "sf")) {
    
    # If not sf object, do these columns exist? If not, stop
    if (!all(c("long", "lat") %in% names(spatial_data))) { 
      stop("Input data must have 'long' and 'lat' columns or be an sf object.")
    }
    sf_data <- sf::st_as_sf(spatial_data, coords = c("long", "lat"), crs = crs)
  } else {
    sf_data <- spatial_data
  }
  
  # Extract coordinates for duplicate checking if it's already sf
  coords <- sf::st_coordinates(sf_data)
  spatial_data$long <- coords[,1]
  spatial_data$lat <- coords[,2]
  
  # Check feature column 
  if (!is.null(feature_col)) {
    if (!(feature_col %in% names(spatial_data))) {
      stop(paste0("Column '", feature_col, "' not found in the input data."))
    }
    
    # Save feature variable
    val_type <- spatial_data[[feature_col]]
    
    # If feature is character or factor, only 1 value per location allowed
    if (is.character(val_type) || is.factor(val_type)) {
      
      # Identify locations with > 1 unique category
      dupes <- aggregate(as.formula(paste(feature_col, "~ long + lat")), 
                         data = spatial_data, 
                         FUN = function(x) length(unique(x)))
      
      if (any(dupes[[feature_col]] > 1)) {
        stop("Categorical duplicates detected: Some GPS coordinates are assigned to multiple different categories. 
             Please resolve these overlaps.")
      }
    }
  }
  
  # Get bounding box and apply buffer
  bbox <- sf::st_bbox(sf_data)
  ext_box <- terra::ext(bbox$xmin - buffer, bbox$xmax + buffer, 
                        bbox$ymin - buffer, bbox$ymax + buffer)
  
  # Create the SpatRaster
  base_rast <- terra::rast(ext_box, 
                           res = res, 
                           crs = crs)
  
  # Optional: populate the raster
  if (!is.null(feature_col)) {
    
    # Check for NA values in the feature column which might skew results
    if (any(is.na(sf_data[[feature_col]]))) {
      message("Note: Feature column contains NA values; these will be ignored in rasterization.")
    }
    
    # Rasterize with features and given function if GPS duplicates 
    out_rast <- terra::rasterize(terra::vect(sf_data), 
                                 base_rast, 
                                 field = feature_col, 
                                 fun = fun)
  } else {
    # Default behavior: Empty raster initialized to 0
    out_rast <- base_rast
    terra::values(out_rast) <- 0
  }
  
  # Handle missing data (uncertainty versus zero)
  if (na_to_zero) {
    message("Empty cells will be turned into zeros.")
    out_rast[is.na(out_rast)] <- 0
  }
  
  return(out_rast)
}
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


