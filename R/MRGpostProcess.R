#' Make some final adjustments to the multiresolution grids
#' 
#' @eval MRGparam("himg")
#' @eval MRGparam("vars")
#' @eval MRGparam("remCols")
#' @eval MRGparam("rounding")
#' 
#' @details
#' The postprocessing function is normally called directly from \code{\link{multiResGrid}}. 
#' However, it might be useful to check the values of the grid cells that
#' will be suppressed, and the values before rounding. In that case 
#' \code{\link{multiResGrid}} can be called with the argument \code{postProcess = FALSE}, 
#' and the post processing be done separately.
#'
#' @examples
#' \donttest{
#' library(sf)
#'
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#'
#' # Set the base resolutions, and create a hierarchical list with gridded data
#' ress = 1000*2^(1:7)
#' ifl = gridData(ifg, "UAA", res = ress)
#' himg = multiResGrid(ifl, ifg = ifg, var = "UAA", weight = "EXT_CORE", postProcess = FALSE)
#' himgp = MRGpostProcess(himg, var = "UAA")
#' 
#' # Confidential grid cells, being suppressed in postProcessing
#' himg[himg$confidential,]
#'  }
#' @export
MRGpostProcess = function(himg, vars, remCols = TRUE, rounding = -1) {
  if (missing(vars) & !is.null(attr(himg, "vars"))) vars = attr(himg, "vars")
  if (missing(remCols) & !is.null(attr(himg, "remCols"))) remCols = attr(himg, "remCols")
  if (missing(rounding) & !is.null(attr(himg, "rounding"))) rounding = attr(himg, "rounding")
  himg[himg$confidential, c("count", "countw")] = NA
  if (!missing(vars) & !isFALSE(rounding)) {
    for (ivar in 1:length(vars)) {
      var = vars[ivar]
      himg[[var]][himg$confidential] = NA
      himg[[paste0("weight", ivar)]][himg$confidential] = NA
      himg[[var]] = round(himg[[var]], rounding)
      himg[[paste0("weight", ivar)]] = round(himg[[paste0("weight", ivar)]], rounding)
    }
  }
  #' @importFrom tidyselect matches
  matchrem = c("small", "reliability", "dom", "freq", "idcount", "idfail", "vres", "idRem", "confidential", "ufun" )
  if (!missing(vars) && length(vars) > 0) { 
    mm = unlist(lapply(1:length(matchrem), FUN = function(x) {ii = grep(matchrem[x], vars); if (length(ii) > 0) x else NULL}))
  } else mm = NULL
  if (length(mm) > 0) matchrem = matchrem[-mm]  
  matchrem = paste(matchrem, collapse = "|")
  matchrem = "small|reliability|idcount|idfail|vres|idRem|confidential|ufun|dom|freq"
  
  if (remCols) himg = himg %>% select(!matches(matchrem))
himg
}
