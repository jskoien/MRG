#' Convenience function based on ggplot2 to plot multi-resolution grids with some
#' default suggestions For full flexibility it is better to use 
#' ggplot2 directly.The function can also be used for ordinary grids 
#' 
#' @eval MRGparam("himg")
#' @eval MRGparam("var") 
#' @eval MRGparam("linecolor")
#' @eval MRGparam("lwd") 
#' @eval MRGparam("borders") 
#' @eval MRGparam("name") 
#' @eval MRGparam("title") 
#' @eval MRGparam("xlim") 
#' @eval MRGparam("ylim") 
#' @eval MRGparam("crs") 
#' @eval MRGparam("clip")  
#' @eval MRGparam("transform") 
#' @eval MRGparam("show.legend") 
#' @eval MRGparam("option") 
#' @eval MRGparam("lwdb") 
#' 
#' 
#' 
#' @details
#' The function is a wrapper around ggplot, possibly calling \code{geom_sf} twice,
#' for the grid itself and for the borders. The function uses the 
#' \code{\link[viridis]{scale_color_viridis}} color scale.
#'   
#' The function will plot the object, and also return a valid ggplot-object that 
#' can be further customized.
#' 
#' @examples
#' 
#' \donttest{
#' library(sf)
#' 
#' if (require(giscoR)) {
#'   useBorder = TRUE 
#' } else {
#'   useBorder = FALSE
#'   print("You need to install giscoR for plotting borders and clipping the gridded maps")
#' }
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' 
#' if (useBorder) {
#' # Read country borders, only used for plotting
#'   borders = gisco_get_nuts(nuts_level = 0)
#' }
#' 
#' ress = c(1,5,10,20,40, 80, 160)*1000
#' # Gridding Utilized agricultural area (UAA)
#' ifl = gridData(ifg, "UAA",res = ress)
#' 
#' # Create a multi-resolution grid of UAA
#' himg1 = multiResGrid(ifl, vars = "UAA", ifg = ifg)
#' 
#' if (useBorder) {
#'   p1 = MRGplot(himg1, UAA, transform = "log10", borders = borders)
#' } else {
#'   p1 = MRGplot(himg1, UAA, transform = "log10")
#' }
#' p1
#' 
#' # Plot can be customized further (reverting to ggplot default color scale in this case)
#' p1 + scale_color_continuous() + scale_fill_continuous()
#' 
#' }
#' 
#' 
#' @export
MRGplot = function(himg, var, linecolor, option = "D", lwd = 0, lwdb = 1, borders, name = waiver(), 
                       title = NULL, xlim, ylim, crs, clip = TRUE, 
                   transform = "identity", show.legend = TRUE) {
  #' @importFrom ggplot2 waiver geom_sf ggplot aes ggtitle theme_bw coord_sf
  if (missing(crs)) crs = st_crs(himg) else himg = st_transform(himg, crs = crs)
  if (missing(var)) {
    stop("A variable name is necessary for plotting")
  } else if (!rlang::quo_is_symbol(rlang::enquo(var))) {
    stop("The variable name must be unqouted - in typical dplyr-style")
  }
  if (!missing(borders)) {
    if (st_crs(himg) != st_crs(borders)) borders = st_transform(borders, crs = st_crs(himg))
    if (clip) himg = st_intersection(himg, borders)
    himg$area = st_area(himg)
    units(himg$area) = NULL
    himg = himg[himg$area > 0, ]
  }
  if (missing(xlim)) xlim = st_bbox(himg)[c(1, 3)]
  if (missing(ylim)) ylim = st_bbox(himg)[c(2, 4)]
  #' @importFrom grDevices colors
  if (!missing(linecolor) && linecolor %in% colors()) {
    p1 = ggplot() + geom_sf(data = himg, aes(fill = {{var}}, color = linecolor), lwd = lwd, show.legend = show.legend) 
  } else if (!missing(linecolor)) {
    p1 = ggplot() + geom_sf(data = himg, aes(fill = {{var}}, color = {{linecolor}}), lwd = lwd, show.legend = show.legend) 
  } else {
    p1 = ggplot() + geom_sf(data = himg, aes(fill = {{var}}, color = {{var}}), lwd = lwd, show.legend = show.legend) 
  }
  #' @importFrom viridis scale_fill_viridis scale_color_viridis
  p1 = p1 + scale_fill_viridis(name = name, trans = transform, option = option) +
  scale_color_viridis(name = name, trans = transform, option = option) +
  ggtitle(title) +
  theme_bw()
  if (!missing(borders)) p1 = p1 + geom_sf(data = borders, fill = NA, colour='black', lwd = lwdb) 
  p1 = p1 + coord_sf(crs = crs, xlim = xlim, ylim = ylim) 
  p1
}


