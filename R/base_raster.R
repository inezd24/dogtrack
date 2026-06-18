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
#'   feature_col = NULL, # Withot a feature
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
      stop(
        "Input data must have 'long' and 'lat' columns or be an sf object."
        )
    }
    
    # Check that there are no NAs present
    if (any(is.na(spatial-data))) {
      stop(
        "The Spatial dataframe contains NA values. Please fill them. "
      )
    }
    
    # Turn into an sf object
    sf_data <- sf::st_as_sf(spatial_data, 
                            coords = c("long", "lat"), 
                            crs = crs
                            )
  } else {
    # Object is already an sf object
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