#' Function to create a gridded (usually multi-resolution grid) from a data.frame or 
#' csv file with information about the corners and resolution, as typically can be 
#' downloaded from Eurostat. The function can also save the grid as a geo-object.
#' 
#' 
#' @eval MRGparam("df")
#' @eval MRGparam("coords")
#' @eval MRGparam("crs")
#' @eval MRGparam("res")
#' @eval MRGparam("coordscale")
#' @eval MRGparam("dsn")
#' @eval MRGparam("driver")
#' @eval MRGparam("layer")
#' @eval MRGparam("Estat")
#' @eval MRGparam("cignore")
#' @eval MRGparam("ellipsisc")
#' 
#' @details
#' This function is mainly for handling csv files downloaded from Eurostat,
#' but can also be used for data from other sources, which adapt the 
#' same csv-convention as Eurostat. 
#' 
#' The Eurostat-files have x- and y-coordinates that have been projected 
#' in the epsg:3035 projection. However, the coordinates show kilometers,
#' not meters, so they have to be multiplied with 1000. Similar data sets
#' might also be offered by other providers. The multiplication can be 
#' done with coordscale, or with \code{Estat = TRUE} (which also sets \code{crs = 3035})
#' 
#' The function will also check the coordinates, if the range of both 
#' x- and y-coordinates are between 360 and 20000, it would often indicate
#' that the coordinates should be multiplied. The function will suggest to 
#' correct this. If the coordinates are actually correct, the 
#' check can be overrun with \code{cignore = TRUE}
#' 
#' If plotting to file, the 
#' 
#' @export
MRGfromDF = function(df, coords = c("x", "y"), coordscale, crs = NA, res = "res", Estat = TRUE, 
                     cignore = FALSE, dsn, layer, ...) {
  #' @importFrom utils read.csv
  if (is.character(df)) df = read.csv(df, ...) else if (!inherits(df, "data.frame")) stop("input must be csv file or data.frame")
  if (Estat & missing(coordscale)) coordscale = 1000 else if (missing(coordscale)) coordscale = 1
  if (Estat & is.na(crs)) crs = 3035
  
  x = coords[1]
  y = coords[2]

  if (coordscale == 1) {
    drx = diff(range(df[[x]]))
    dry = diff(range(df[[y]]))
    
    if (drx > 360 & drx < 20000 & dry > 360 & dry < 20000 &!cignore) {
      cat("The values of the coordinates seem neither to be lat-lon or projected\n
          Could it be that the coordinates should be scaled (argument coordscales)?\n")
      response <- readline("Give coordscale or answer NO")
      if (is.numeric(response))  coordscale = response
    }
  }  
  if (coordscale != 1) df = df %>% mutate(!!x := !! as.name(x)*coordscale, !!y := !! as.name(y)*coordscale)
  
  # Function to create a square polygon
  create_square <- function(x, y, res) {
    # Define the coordinates of the square
    coords <- matrix(c(x, y, 
                       x + res, y, 
                       x + res, y + res, 
                       x, y + res, 
                       x, y), ncol = 2, byrow = TRUE)
    
    # Create a polygon from the coordinates
    #' @importFrom sf st_polygon st_sf st_crs
    polygon <- st_polygon(list(coords))
    return(polygon)
  }
  
  # Apply the function to your data.frame
#  df$polygons <- mapply(create_square, df[[x]], df[[y]], df[[res]], simplify = FALSE)
#  df$polygons <- lapply(1:dim(df)[1], FUN = function(ii) create_square(df[[x]][ii], df[[y]][ii], df[[res]][ii]))
  pols <- lapply(1:dim(df)[1], FUN = function(ii) create_square(df[[x]][ii], df[[y]][ii], df[[res]][ii]))
  #
  # Convert the data.frame to an sf object
  sf_df <- st_sf(df, geometry = pols)
#' @importFrom sf st_crs<- st_write
  if (!missing(crs)) st_crs(sf_df) = st_crs(crs)
  # Remove the x, y, res, and polygons columns
  sf_df = sf_df[, -which(names(sf_df) %in% c(x, y))]
  if (!missing(dsn)) if (missing(layer)) st_write(sf_df, dsn) else st_write(sf_df, dsn, layer)
  invisible(sf_df)
 }

  
  
