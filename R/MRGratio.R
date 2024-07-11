#' This function will compute the ratio between two variables
#' 
#' 
#' @param himg A highest possible resolution grid, created with the multiResGrid function, for a variable
#' @param ifg An sf-object with the locations and the data of the holdings
#' @param var1 The variable of the nominator in the ratio
#' @param var2 The variable of the denominator in the ratio 
#' @param weight1 The weights of the variable of the nominator in the ratio
#' @param weight2 The weights of the variable of the denominator in the ratio 
#' @param mincount The minimum number of farms for a grid cell (threshold rule)
#'
#' @param nlarge Parameter to be used if the nlarge(st) farms should count for maximum plim percent of
#'         the total value for the variable in the grid cell (see details)
#' @param plim See nlarge, and details below
#' @param domEstat Should the dominance rule be applied as in the IFS handbook (TRUE), where 
#'          the weights are rounded before finding the first nlarge contributors, or should 
#'         it be the first nlarge contributors*weight, where also fractions are considered (FALSE)? 
#' @eval MRGparam("postProcess")
#' @param rounding rounding is usually not applied for ratios, hence the high default value (meaning no rounding)"
#' 
#' @examples
#' library(sf)
#' library(viridis)
#' library(ggplot2)
#' library(patchwork)
#' library(giscoR)
#'
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' # Read country borders, only used for plotting
#' borders = gisco_get_nuts(nuts_level = 0)
#' dkb = borders[borders$CNTR_CODE == "DK",] %>% st_transform(crs = 3035)
#' 
#' ifg_org = ifg[ifg$UAAXK0000_ORG > 0, ]
#' ress = c(1,5,10,20,40,80)*1000
#' ifl = gridData(ifg_org, "UAAXK0000_ORG", res = ress)
#' # Create multi-resolution grid of organic UAA
#' himg1 = multiResGrid(ifl, ifg = ifg_org, var = "UAAXK0000_ORG")
#' 
#' # Create the ratio between organic UAA and total UAA per grid cell
#' himg2 = MRGratio(himg1, ifg, "UAAXK0000_ORG", "UAA")
#' g2 = ggplot() + geom_sf(data = himg2, aes(fill = ratio*100)) +
#'   scale_fill_viridis( name = bquote(atop("Percentage Org UAA", "of total UAA"))) +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle(bquote("Ratio between organic UAA and total UAA")) + 
#'   theme_bw()
#' g2
#' 
#' 
#' # Create the ratio between organic UAA and total UAA per grid cell
#' himg2 = MRGratio(himg1, ifg, "UAAXK0000_ORG", "UAA")
#' g2 = ggplot() + geom_sf(data = himg2, aes(fill = ratio*100)) +
#'   scale_fill_viridis( name = bquote(atop("Percentage Org UAA", "of total UAA"))) +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle(bquote("Ratio between organic UAA and total UAA")) + 
#'   theme_bw()
#' g2
#' 
#' 
#' @export
MRGratio = function(himg, ifg, var1, var2, weight1, weight2 = weight1, mincount = 10, nlarge = 2, plim = 0.85,
                    domEstat = TRUE, postProcess = TRUE, rounding = 10) {
  if (missing(weight1)) ifg$w1 = 1; weight1 = weight2 = "w1"
  himgdat = sf::st_drop_geometry(himg)
  if (!"ID" %in% names(ifg)) ifg$ID = 1:dim(ifg)[1]
  loh = as.data.frame(st_join(ifg, himg, join = st_within))
  names(loh) = gsub(".x", "", names(loh), fixed = TRUE)
  #' @importFrom stats aggregate
  agr = aggregate(loh[, var2]*loh[, weight2], by = list(loh$ID.y), FUN = sum)
  himg$newvar = agr$x[match(himg$ID, agr$Group.1)]  
  himg$ratio = himgdat[,var1]/himg$newvar
  himg$confidentiality = FALSE
# Check dominance for newvar
  for (ID in himg$ID) {
    hIDs = loh$ID[which(loh$ID.y == ID)]
    ifglldat = loh[loh$ID %in% hIDs,]
    ifglldat$weight = ifglldat[, weight2]
    ifglldat$gridvar = ifglldat[, var2]
    hcs = sum(ifglldat[,weight2])
    if (hcs < mincount) {
      himg$confidentiality[ID] = TRUE
    } else {
      himg$confidentiality[ID] = dominanceRule(ifglldat, nlarge, plim, domEstat = domEstat)
    }      
  }
  if (postProcess) himg = MRGpostProcess(himg, "ratio", remCols = FALSE, rounding = rounding)
  himg
}
