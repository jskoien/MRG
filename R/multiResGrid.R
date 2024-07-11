#' Function that creates a multi-resolution grid with larger grid cells in
#' regions with lower resolution of farms, or where farms need to
#' be anonymized for disclosure control reasons
#'
#' Two main confidentiality rules are considered:
#' - Threshold rule (suppression due to a minimum number of counts)
#' - Dominance rule (suppression due to dominance by one or more units)
#'
#' @eval MRGparam("MRGobject")
#' @eval MRGparam("gdl")
#' #@eval MRGparam("ress")
#' @eval MRGparam("mincount")
#' @eval MRGparam("nlarge")
#' @eval MRGparam("plim")
#' @eval MRGparam("ifg")
#' @eval MRGparam("vars")
#' @eval MRGparam("weights")
#' @eval MRGparam("countFeatureOrTotal")
#' # @eval MRGparam("minpos")
#' @eval MRGparam("verbose")
#' @eval MRGparam("nclus")
#' @eval MRGparam("clusType")
#' @eval MRGparam("plim")
#' @eval MRGparam("domEstat")
#' @eval MRGparam("outfile")
#' @eval MRGparam("checkDominance")
#' @eval MRGparam("checkReliability")
#' @eval MRGparam("userfun")
#' @eval MRGparam("strat")
#' @eval MRGparam("confrules")
#' @eval MRGparam("suppresslim")
#' @eval MRGparam("sumsmall")
#' @eval MRGparam("suppresslimSum")
#' @eval MRGparam("reliabilitySplit")
#' @eval MRGparam("plotIntermediate")
#' @eval MRGparam("addIntermediate")
#' @eval MRGparam("postProcess")
#' @eval MRGparam("rounding")
#' @eval MRGparam("remCols")
#' @eval MRGparam("ellipsis")
#'
#' @details This function will find the highest resolution data set that fulfills the 
#'        confidentiality rules and potential reliability rules for variable(s) of interest.
#'        Starting with the second highest resolution (5 km in the default settings),
#'        the function will check if any of the 1 km sub pixels will have values not fulfilling
#'        any of the confidentiality rules (number of farms, values of the 2 largest compared to values of the 
#'        entire grid cell).
#'        If all values are above the confidentiality limits, the grid cells will be kept at a 1 km resolution,
#'        otherwise only the 5 km grid cell will be kept. This will again be tested against the confidentiality
#'        rules in the next iteration, when grid cells will possibly be merged to 10 km grid cells.
#'
#'        The standard threshold rule for spatial data is at least 10 units (mincount). 
#'        The parameters nlarge and plim are used for determining the dominance treatment for the variable of interest,
#'        with default values of nlarge=2 and plim=0.85. 
#'        If more than plim of the values of the grid cell (e.g. UAA, arable land, number of livestock)
#'        is explained by 1-nlarge weighted holdings, the grid cell will not pass the confidentiality rule.
#'        
#'        The concept of reliability is explained in details in section 4.6 in the integrated farm survey handbook for 2023:
#'        https://wikis.ec.europa.eu/display/IFS/Integrated+Farm+Statistics+Manual+%7C+2023+edition
#'        In short, it is an estimate of the coefficient of variation for an estimate (a grid cell in this case),
#'        based on the number in the sample relative to the number in the population, and taking into account
#'        possible stratified sampling approaches. The number is zero if all holdings in the population in 
#'        a grid cell has been sampled, and the default requirement is that the CV is less than 35%.
#'        
#'        There are some cases where aggregation might not be desired. In the situation where a 
#'        relatively large single grid cell does not respect the confidentiality rules, it is fine to 
#'        aggregate it if the neighbouring grid cells are also relatively large. However, it can be seen
#'        as unfortunate if the single cell was aggregated with many smaller grid cells that could otherwise 
#'        be disseminated at a high resolution. The added value of being able to present a value for a 
#'        region with very few farms is perhaps lower than what is lost by having to aggregate to a 
#'        lower resolution. The parameter \code{suppresslim} indicates the minimum value in a grid cell 
#'        relative to the possible lower resolution grid cell 
#'        before it is necessary to aggregate. If the limit is 0.05, a grid cell would only cause an aggregation
#'        to lower resolution if the value in the grid cell is more than 5% of the value in the lower resolution 
#'        grid cell. Instead, it would be left as it is, and will be suppressed in the post-processing step.
#'
#'        There are cases when the built-in confidentiality checks are not what the user needs. 
#'        That is why it is possible to submit a user defined function. This function needs to follow
#'        certain rules. 
#'        \enumerate{
#'          \item The first argument must be a data.frame with name \code{df}. 
#'                This is a data.frame with the individual records
#'                for a particular grid cell. It has three columns: 
#'              \enumerate{
#'                 \item himgid - the ID of the current grid cell. This is the grouping variable
#'                        and is constant for the data.frame
#'                 \item gridvar - a new common name for the current variable to be gridded
#'                 \item weight - the weight of the variable to be gridded
#'              }  
#'          \item The function can include additional parameters for calculation of confidentiality
#'                 (or reliability, or suitability, if the meaning of the function refers to something else).
#'                 This can be new parameters to this particular function (through the 
#'                 ellipsis argument (...) of \code{multiResGrid}), existing parameters to \code{multiResGrid},
#'                 or potentially internal variables of \code{multiResGrid.})
#'          \item The result of the function must be a logical, either the rule was passed 
#'                for the records of this grid cell, or not (TRUE/FALSE)
#'          \item The function can potentially use internal variables of \code{multiResGrid}, however,
#'                the meaning of these will have to be understood from the code
#'        }
#'        A simple example of a \code{userfun} is given in the example section below (the one producing \code{himg6})
#'        
#'
#' @examples
#' \donttest{
#' library(sf)
#' library(viridis)
#' library(ggplot2)
#' library(patchwork)
#' library(giscoR)
#'
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' ifs_weight = ifs_dk %>% dplyr::filter(Sample == 1) # Extract weighted subsample
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' fsg = fssgeo(ifs_weight, locAdj = "LL")
#' # Read country borders, only used for plotting
#' borders = gisco_get_nuts(nuts_level = 0)
#' dkb = borders[borders$CNTR_CODE == "DK",] %>% st_transform(crs = 3035)
#'
#' # Set the base resolutions, and create a hierarchical list with gridded data
#' ress = c(1,5,10,20,40, 80, 160)*1000
#' # Gridding Utilized agricultural area (UAA)
#' ifl = gridData(ifg, "UAA",res = ress)
#' # Gridding organic utilized agricultural area
#' ifl2 = gridData(ifg, vars = "UAAXK0000_ORG", res = ress)
#' 
#' # Gridding UAA and organic UAA together
#' ifl3 = gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), res = ress)
#' 
#' # Gridding the UAA from the survey - the survey weights are in the column EXT_MODULE
#' fsl = gridData(fsg,  vars = c("UAA"), weights = "EXT_MODULE",  res = ress)
#' 
#' # Create a multi-resolution grid only with farm number as confidentiality rule, then plot results
#' himg0 = multiResGrid(ifl, checkReliability = FALSE, suppresslim = 0)
#' ggplot(himg0) + geom_sf(aes(fill = count))
#' 
#' # Create a multi-resolution grid of UAA, also based on the dominance rule (default)
#' himg1 = multiResGrid(ifl, vars = "UAA", ifg = ifg)
#' ggplot(himg1) + geom_sf(aes(fill = UAA))
#' 
#' # Create multi-resolution grid of organic UAA
#' himg2 = multiResGrid(ifl2, vars = "UAAXK0000_ORG", ifg = ifg)
#' himg21 = multiResGrid(ifl2, vars = "UAAXK0000_ORG", ifg = ifg, postProcess = FALSE)
#' 
#' ggplot(himg2) + geom_sf(aes(fill = UAAXK0000_ORG))
#' 
#' # Create joint multi-resolution grid of organic UAA and total UAA
#' himg3 = multiResGrid(ifl3, vars = c("UAA", "UAAXK0000_ORG"), ifg = ifg, 
#'                   checkReliability = FALSE, suppresslim = 0)
#' p1 = ggplot(himg3) + geom_sf(aes(fill = UAA))
#' p2 = ggplot(himg3) + geom_sf(aes(fill = UAAXK0000_ORG))
#' p1 + p2
#' 
#' # Create multi-resolution grid of UAA, based on survey data,
#' # with and without applying reliability check
#' # Slow!
#' himg4 = multiResGrid(fsl,  vars = c("UAA"), weights = "EXT_MODULE", ifg = fsg, 
#'                       strat = "STRA_ID_CORE", checkReliability = FALSE)
#' himg5 = multiResGrid(fsl,  vars = c("UAA"), weights = "EXT_MODULE", ifg = fsg, 
#'                       strat = "STRA_ID_CORE", checkReliability = TRUE)
#'                       
#'# Apply suppreslim to suppress insignificant grid cells
#'# Show intermediate maps of confidential cells (wait 5 seconds)
#' pint = ifelse(interactive(), 5, FALSE)
#' himg11 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0, plotIntermediate = pint)
#' himg12 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0.02, plotIntermediate = pint)
#' himg13 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0.05, plotIntermediate = pint)
#' himg14 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0.1, plotIntermediate = pint)
#' himg15 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0.2, plotIntermediate = pint)
#'  
#'  
#'  # This is an example of a userfun that can be used for alternative restrictions
#'  # for a grid cell. This particular toy example assures that there are at least
#'  # nabove records with a value (UAA in this case) above limit. 
#'  ufun = function(df, nabove, limit) {
#'    sum(df$gridvar > limit) < nabove
#'  }
#'  
#' himg6 = multiResGrid(ifl, vars = "UAA", ifg = ifg, 
#'                  suppresslim = 0.2, plotIntermediate = pint, userfun = ufun, nabove = 5, limit = 10)
#'  
#'  
#'  
#' himg00 = st_intersection(dkb, himg0)
#' ggplot() + geom_sf(data = himg00, aes(fill = count, color = count)) +
#'   scale_fill_viridis( name = "number of farms", trans = "log10") +
#'   scale_color_viridis( name = "number of farms", trans = "log10") +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("Number of farms for variable grid cell size, only frequency confidentiality") +
#'   theme_bw()
#'  
#' himg01 = st_intersection(dkb, himg1)
#' ggplot() + geom_sf(data = himg01, aes(fill = count, color = count)) +
#'   scale_fill_viridis( name = "number of farms", trans = "log10") +
#'   scale_color_viridis( name = "number of farms", trans = "log10") +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("Number of farms for variable grid cell size, frequency and dominance confidentiality") +
#'   theme_bw()
#' 
#'   
#' # Plot the density of organic agriculture, as hectares per square km
#' himg02 = st_intersection(dkb, himg2)
#' himg02$orgarea = himg02$UAAXK0000_ORG/units::set_units(st_area(himg02), "km^2")
#' units(himg02$orgarea) = NULL
#' ggplot() + geom_sf(data = himg02, aes(fill = orgarea), lwd = 0) +
#'   scale_fill_viridis( name = "ha / km2") +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("Organic UAA density")  +
#'   theme_bw()
#' 
#' # Plot the relative abundance of organic UAA relative to total UAA
#' himg03 = st_intersection(dkb, himg3)
#' himg03$ouaashare = himg03$UAAXK0000_ORG/himg03$UAA*100
#' ggplot() + geom_sf(data = himg03, aes(fill = ouaashare), lwd = 0) +
#'   scale_fill_viridis( name = "% Organic") +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("Organic share")  +
#'   theme_bw()
#'   
#'   
#' # Plot maps from survey data before and after adding the reliability constraint 
#' # The percentage of UAA can be above 100% due to farm area being registered at the location
#' # of the administration building, but the map without reliability check has too high values 
#' # for too many cells
#' 
#' himg04 = st_intersection(dkb, himg4)
#' himg04$area = st_area(himg04)/1e6
#' units(himg04$area) = NULL
#' himg04$uaashare = himg04$UAA/himg04$area
#' himg04$uaashare[himg04$uaashare > 1000] = 1000
#' g4 = ggplot() + geom_sf(data = himg04, aes(fill = uaashare), lwd = 0) +
#'   scale_fill_viridis( name = "% UAA",  trans = "log10", limits = c(1,1000)) +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("UAA share (sample without reliability check)")  +
#'   theme_bw()
#'   
#' himg05 = st_intersection(dkb, himg5)
#' himg05$area = st_area(himg05)/1e6
#' units(himg05$area) = NULL
#' himg05$uaashare = himg05$UAA/himg05$area
#' himg05$uaashare[himg05$uaashare > 1000] = 1000
#' g5 = ggplot() + geom_sf(data = himg05, aes(fill = uaashare), lwd = 0) +
#'   scale_fill_viridis( name = "% UAA",  trans = "log10", limits = c(1,1000)) +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("UAA share (sample with reliability check)")  +
#'   theme_bw()
#'   
#'g4 + g5 + plot_layout(guides = "collect")
#'   
#' himg06 = st_intersection(dkb, himg6)
#' ggplot() + geom_sf(data = himg06, aes(fill = UAA), lwd = 0) +
#'   scale_fill_viridis( name = "ha") +
#'   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
#'   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
#'   ggtitle("UAA, with additional user defined function")  +
#'   theme_bw()
#'   
#'      
#' # Plot the different maps from using different suppreslim values
#' himgs = list(himg11, himg12, himg13, himg14, himg15)
#' slims = c(0, 0.02, 0.05, 0.1, 0.2)
#' plots = list()
#' uaas = c(himg11$UAA, himg12$UAA, himg13$UAA, himg14$UAA, himg15$UAA)
#' lims = range(uaas[uaas > 0], na.rm = TRUE)
#' for (ii in 1:5) {
#'   himg = st_intersection(dkb, himgs[[ii]])
#'   plots[[ii]] = 
#'    ggplot() + geom_sf(data = himg, aes(fill = UAA), lwd = 0) +
#'     scale_fill_viridis( name = "UAA (ha)", trans = "log10", limits = lims, na.value="red") +
#'     geom_sf(data = dkb, fill = NA, colour='black', lwd = 0.5) +
#'     ggtitle(paste("Suppresslim = ", slims[[ii]])) +
#'     xlab("") + ylab("") +
#'     theme_bw()
#' }
#' 
#' plots[[1]]  + plots[[2]] + plots[[3]]  + plots[[4]] + plots[[5]] + plot_layout(guides = "collect")
#'  
#' 
#' }
#' 
#' 
#' @rdname multiResGrid
#' @export
multiResGrid.MRG <- function(MRGobject, ...) {
  dots = list(...)
  #' @importFrom utils modifyList
  if (length(dots) > 0) MRGobject = modifyList(MRGobject, dots)
  do.call(multiResGrid, MRGobject)
}
#'
#' @rdname multiResGrid
#' @export
multiResGrid.list <- function(gdl, ifg, vars, weights, countFeatureOrTotal = "feature", mincount = 10, #minpos = 4, 
                       nlarge = 2,
                          plim = 0.85, verbose = FALSE, nclus = 1, clusType, domEstat = TRUE, 
                          outfile = NULL, checkDominance = TRUE,
                          checkReliability = FALSE, userfun, strat = NULL, confrules = "individual", 
                          suppresslim = 0, sumsmall = FALSE, suppresslimSum = NULL,
                          reliabilitySplit = TRUE,
                          plotIntermediate = FALSE,  addIntermediate = FALSE,
                          postProcess = TRUE, rounding = -1, remCols = TRUE, ...) {
#  To avoid R CMD check notes
  hsum = small = weight = data = himgid = dominance = . = NULL
  if (!missing(ifg) && !inherits(ifg, "sf")) stop("ifg is not an sf-object ")
  ress = unlist(lapply(gdl, FUN = function(gdll) gdll$res[1]))
  
  if (checkReliability) {
    if (missing(strat) | is.null(strat) ) {
      if (!"strat" %in% names(ifg)) ifg$strat = 1  
      #' @importFrom dplyr mutate group_by
    } else ifg = ifg %>% mutate(strat = .data[[strat]])
  }
  if (!missing(vars)) {
    if (missing(ifg))
      stop(paste("Cannot create values for variable(s) ",
                 vars, " without ifg"))
    
    # Some kind of test needed to check if this step has already been done, and also if it is necessary    
    #    for (iw in 1:length(vars)) ifg[, paste0(vars, iw)] = data.frame(ifg)[, vars[iw]]
    if (!"ID" %in% names(ifg)) ifg$ID = 1:dim(ifg)[1]
    for (iw in 1:length(vars)) ifg[, paste0("gridvar", iw)] = st_drop_geometry(ifg[, vars[iw]])
    ifg = addweights(ifg, vars, weights)
    wts = paste0("weight", 1:length(vars))
    if (checkReliability) {
      ifg = ifg[, c("ID", paste0("gridvar", 1:length(vars)), paste0("weight", 1:length(vars)), "strat")]
    } else {
      ifg = ifg[, c("ID", paste0("gridvar", 1:length(vars)), paste0("weight", 1:length(vars)))]
    }
  } 
  
  himg = gdl[[1]] %>% mutate(confidential = FALSE, reliability = FALSE, small = FALSE,
                               freq = FALSE, dom = FALSE, ufun = FALSE)
  
  himgs = list()
  lohs = list()
  if (!missing(vars)) for (ivar in 1:length(vars)) himg[,paste0("vres", ivar)] = 0
  for (ires in 2:(length(ress) + 1)) {
    lres = ress[ires]
    if (ires <= length(ress)) {
      limg = gdl[[ires]] %>% mutate(confidential = FALSE, reliability = FALSE, small = FALSE,
                                      freq = FALSE, dom = FALSE, ufun = FALSE)
      if (!missing(vars)) for (ivar in 1:length(vars)) limg[,paste0("vres", ivar)] = 0
    }    
    ##    #' @importFrom utils txtProgressBar setTxtProgressBar
    #' @importFrom sf st_join st_within st_drop_geometry st_make_grid st_coordinates st_crs st_area
    if (!missing(ifg) && !is.null(ifg) && ires <= length(ress)) {
      ifg$himgid = st_join(ifg, himg, join = st_within)$ID.y
      ifg$limgid = st_join(ifg, limg, join = st_within)$ID.y
      ifgl = data.frame(ifg)[, grep("ID|himgid|gridvar|weight|strat", names(ifg))]
    } else {
      ifgl = NULL
      if (!missing(ifg)) {
        ifg$himgid = st_join(ifg, himg, join = st_within)$ID.y
        ifgl = data.frame(ifg)[, grep("ID|himgid|gridvar|weight|strat", names(ifg))]
      } 
    }
    if (ires <= length(ress)) loh = st_drop_geometry(st_join(himg, limg, join = st_within))
    if (!missing(vars) & ires <= length(ress)) {
      for (ivar in 1:length(vars)){
        if (tolower(countFeatureOrTotal) == "feature") {
#          sel = (loh[[paste0(vars[ivar], "_w", ivar, ".x")]]  <
#                   suppresslim*loh[[paste0(vars[ivar], "_w", ivar, ".y")]]) & loh[[paste0("weight",ivar, ".x")]] < mincount
          sel = (loh[[paste0(vars[ivar], ".x")]]  <
                 suppresslim*loh[[paste0(vars[ivar], ".y")]]) & loh[[paste0("weight",ivar, ".x")]] < mincount
#        } else sel = (loh[[paste0(vars[ivar], "_w", ivar, ".x")]]  <
#                        suppresslim*loh[[paste0(vars[ivar], "_w", ivar, ".y")]]) & loh[["countw.x"]] < mincount
      } else sel = (loh[[paste0(vars[ivar], ".x")]]  <
                        suppresslim*loh[[paste0(vars[ivar], ".y")]]) & loh[["countw.x"]] < mincount
        if (sumsmall & sum(sel) == 0) {
          sshare = data.frame(Group.1 = 999999, x = 9999999)
        } else if (sumsmall) {
#          sshare = aggregate(loh[[paste0(vars[ivar], "_w", ivar, ".x")]][sel], by = list(loh$ID.y[sel]), sum) 
          sshare = aggregate(loh[[paste0(vars[ivar], ".x")]][sel], by = list(loh$ID.y[sel]), sum) 
          loh$sshare = sshare$x[match(loh$ID.y, sshare$Group.1)]
#          loh$sshare[is.na(loh$sshare)] = max(loh[[paste0(vars[ivar], "_w", ivar, ".y")]])
          loh$sshare[is.na(loh$sshare)] = max(loh[[paste0(vars[ivar], ".y")]])
#          loh[[paste("hsmall", ivar)]] = loh$sshare <  suppresslimSum*loh[[paste0(vars[ivar], "_w", ivar, ".y")]]
          loh[[paste("hsmall", ivar)]] = loh$sshare <  suppresslimSum*loh[[paste0(vars[ivar], ".y")]]
        } else {
          loh[[paste("hsmall", ivar)]] = sel
        }
      }
      
      #' @importFrom tidyselect contains
      himg$small = loh %>% select(contains("hsmall")) %>% mutate(hsum = rowSums(.)) %>%
        mutate(small = (hsum == length(vars))) %>% select(small) %>% pull
    }
    himgdat = st_drop_geometry(himg)
    ifgdat = st_drop_geometry(ifgl)
    
    if (tolower(countFeatureOrTotal) == "feature" & !missing(vars)) {
      ww = himgdat[,names(himgdat) %in% paste0("weight", 1:length(vars)), drop = FALSE]
      ww = apply(ww, MARGIN = 1, FUN = function(x) if (sum(x > 0))  min(x[x>0]) else 0)
    } else ww = himgdat[, "countw"]
    wf = which(ww > 0 & ww < mincount)
    if (length(wf) > 0)  himg$freq[which(himg$ID %in% himgdat$ID[wf])] = TRUE
    len = ifelse(missing(vars), 1, length(vars))
    ehimgid = ufres = NULL
    for (ivar in 1:len) {
      ifgdatl = NULL
      if (checkDominance & !missing(vars)) {
        ifgdatl <- ifgdat[,c("himgid", paste0("gridvar", ivar), paste0("weight",ivar))] 
        names(ifgdatl) = c("himgid", "gridvar", "weight")
        ifgdatl$ehimgid = ifgdatl$himgid
        ifgdatl = ifgdatl[order(ifgdatl$ehimgid),]
        #' @importFrom tidyr unnest nest 
        #' @importFrom purrr map 
        dom = ifgdatl %>% filter(weight != 0)  %>%
          group_by(ehimgid) %>%  nest() %>%
          mutate(dominance = map(data, ~dominanceRule(., nlarge = nlarge, plim = plim, 
                                                       domEstat = domEstat))) %>%
          unnest(dominance) %>% ungroup %>% select(dominance) %>% pull 

        domid = which(dom)
        if (length(domid) > 0) himg$dom[domid] = TRUE
      } 
      if (!missing(userfun) && is.function(userfun)) {
        if (is.null(ifgdatl)) ifgdatl <- ifgdat[,c("himgid", paste0("gridvar", ivar), 
                                                   paste0("weight",ivar))] 
        
        ifgdatl$ehimgid = ifgdatl$himgid
        ifgdatl = ifgdatl[order(ifgdatl$ehimgid),]
        dots = list(...)
          fargs = names(formals(userfun))
          if ("hareas" %in% fargs) hareas = st_area(himg)
          passed_in  <- setdiff(intersect(ls(),fargs),"df")
          names(passed_in) <- passed_in
          args_to_pass <- map(passed_in, dynGet)
          ddots = dots[names(dots) %in% fargs]
          args_to_pass <- c(args_to_pass,ddots)
          
          localUserFun = function(subdata, args_to_pass) {
            local_args <- c(list("df" = subdata), args_to_pass)
            do.call(userfun, args=local_args)}
          
          ufRes = ifgdatl %>% 
            group_by(ehimgid) %>%  nest() %>%
            mutate(ufres = map(data, ~localUserFun(., args_to_pass))) %>%
            unnest(ufres) %>% ungroup %>% select(ufres) %>% pull 
        ufid = which(ufRes)
        if (length(ufid) > 0) himg$ufun = TRUE
      }
    }
    
    if (checkReliability) {
      if (reliabilitySplit & (dim(ifg)[1] > 50000 | dim(himg)[1] > 1000)) {
        if (is.logical(reliabilitySplit)) reliabilitySplit = dim(ifg)[1] %/% 30000
        if (dim(himg)[1] > reliabilitySplit*1000) reliabilitySplit = dim(himg)[1] %/% 1000
      }
      if (!missing(vars) && !is.null(vars)){
        for (ivar in 1:length(vars)){
          vestres = mrg_varestim(ifg, var = paste0("gridvar", ivar), strat = "strat", PSU = "ID", 
                                 weight = paste0("weight", ivar), split = reliabilitySplit)
          himg[,paste0("vres",ivar)] = vestres$rse
        }
      }
      nonvalids = which(apply(st_drop_geometry(himg[, grep("vres", names(himg))]), 1, max, na.rm = TRUE) > 0.35)
      himg$reliability[nonvalids] = TRUE
    }
    
    himg$confidential = rowSums(st_drop_geometry(himg[,c("freq", "dom", "ufun", "reliability")])) > 0
    
    if (ires <= length(ress)) {
      loh$confidential = himg$confidential
      loh$small = himg$small
      idRem = which(himg$confidential & !himg$small)
      # Compute how many non-confidential himg-cells there are per limg-cell
      nclimg = aggregate(!loh$confidential, by = list(loh$ID.y), sum)
      loh$nclimg = nclimg$x[loh$ID.y]
      # Add cells to the aggregation list if none of the hing-cells in a limg-cell
      # are non-confidential, even if tthey are are small
      idAdd = NULL
      idRem = unique(c(idRem, which(loh$small & loh$nclimg == 0)))
      if (length(idRem) > 0) {
        idAdd = which(limg$ID %in% unique(loh$ID.y[loh$ID.x %in% idRem]))
        # Check how many himg-cells per limg-cell. Do not aggregate if it is only 1 
        iac = aggregate(rep(1, length(loh$ID.x)), by = list(loh$ID.y), FUN = sum)
        singlimg = iac$Group.1[iac$x == 1]
        if (length(singlimg) > 0) idAdd = idAdd[!(idAdd %in% singlimg)]
        # Find all himg-cells in the limg-cells to be added 
        idRem = which(himg$ID %in% unique(loh$ID.x[loh$ID.y %in% idAdd]))
        if (length(idRem) == 0) break
        remIdRem = loh$ID.x[loh$ID.y %in% singlimg]
        if (length(remIdRem) > 0) {
          if (sum(idRem %in% remIdRem) > 0) {
            print(sum(idRem %in% remIdRem))
            stop("There are still single himg-cells to be deleted from limg-cells
                This should not happen")
          }
          idRem = idRem[!(idRem %in% remIdRem)]
        }
        
        himg = himg[-idRem,]
        himg = rbind(himg, limg[idAdd,])
      } 
    } else {idRem = NULL; idAdd = NULL}      
    # Reindex himg, and add himg and loh to the lists
    himg$ID = 1:dim(himg)[1]
    himgs[[ires]] = himg
    lohs[[ires]] = loh
    
    print(paste("ires", ires, ress[ires], "#himg-cells:", dim(himg)[1], "; removed:", 
                length(idRem), "; added:", length(idAdd), "; confidential:", sum(himg$confidential) ))
    if (plotIntermediate) {
      plot(himg[, "confidential"], main = "Confidential cells")
      Sys.sleep(ifelse(is.logical(plotIntermediate), 5, plotIntermediate))
    }
  }
  if ("gridvar" %in% names(himg)) himg = himg[, -grep("gridvar", names(himg))]
  if (!missing(vars)) attr(himg, "vars") = vars
  if (postProcess) himg = MRGpostProcess(himg, vars, remCols, rounding)
  if (addIntermediate) {
    attr(himg, "himgs") = himgs
    attr(himg, "lohs") = lohs
  }
  himg
}









mrg_varestim <- function(x, var, strat, PSU, weight, split){
  ID = n = hld = w_sum = NULL
  if (inherits(x, "sf")) x = st_drop_geometry(x)
  if (!missing(PSU)) x$ID = x[[PSU]]  
  if (missing(strat) || is.null(strat)) {
    x$strat = 1
  } else {
    if (!strat %in% names(x)) stop(paste(strat, "is missing from from the data.frame"))
    x$strat = x[[strat]]  
  }
  if (split == 1) {
    df = x
    out_var = vardom(dataset = df, Y= var, H = "strat",
                     PSU = "ID",
                     w_final = weight, Dom = "himgid")$all_result
  } else {
    himgid <- unique(x$himgid)    
    #' @importFrom dplyr left_join mutate ungroup group_by distinct case_when n
    #' @importFrom sjmisc split_var 
    df_cl <- left_join(x, data.frame(himgid = himgid, cluster = split_var(himgid, n = split)))
    out_var <- NULL
    for (isp in 1:split){
      df = x[which(df_cl$cluster == isp),]
      print(paste("df", paste(dim(df)[1], length(unique(df$himgid)))))
      t <- df %>% group_by(ID, strat, himgid) %>% 
        summarise(hld=n(), w_sum = sum(.data[[weight]], na.rm = T)) %>% 
        filter(hld == 1 & w_sum > 1) %>% ungroup
      if (!is.null(t)){
        h_st <- t %>% distinct(ID) %>% pull()
        df <- df %>% mutate(strat = case_when(ID %in% c(h_st)~99999,
                                              T ~ strat))}
      #' @importFrom vardpoor vardom
      est <- vardom(dataset = df, Y= var, H = "strat",
                    PSU = "ID",
                    w_final = weight, Dom = "himgid")
      out_var<-rbind(out_var, est$all_result)
    }
  }
  out_var$himgid = as.numeric(out_var$himgid)
  out_var = out_var[order(out_var$himgid),]
  out_var
}





dominanceRule = function(ifglldat, nlarge, plim, domEstat = TRUE) {
  #' @importFrom dplyr summarise
  Y <- sum(ifglldat$gridvar*ifglldat$weight)
  ifglldat <- ifglldat[order(ifglldat$gridvar, ifglldat$weight, decreasing = TRUE),]
  if (domEstat) {
    # Extrapolated aggregated value of the cell Y
    # Need to loop to account for different nlarge values, even though nlarge=2 is the standard rule
    dominance = FALSE
    nlarge = min(dim(ifglldat)[1], nlarge)
    for (nc in 1:nlarge){
      #' @importFrom magrittr "%>%"
      #' @importFrom dplyr slice select pull arrange ungroup
      #' @importFrom rlang .data
      # wmax= nc largest contributor
      #      wmax<-ifglldat %>% slice(1:nc) %>% select(.data$weight) %>%  pull()
      wmax<-ifglldat$weight[1:nc]
      # This should only round the weights that are above 0.5, to avoid creation of zero weights
      # This should avoid any issues due to low weights, such as
      # wmax = c(0.1, 0.4, 0.7, 1.2, 1.4)
      wmaxr<- ifelse(wmax > 0.5, round(wmax), wmax)
      #      wmaxr<-wmax %>% ifelse(.data > 0.5, round(.data), .data)  #JON: Separated wmax from wmaxr, as only the second should be rounded
      # xmax nc largest values of the variable
      #      xmax<-ifglldat %>% slice(1:nc) %>% select(.data$gridvar) %>% pull()
      xmax<-ifglldat$gridvar[1:nc] 
      # sum of n largest contributor is <= nlarge and aggregated extrapolated value of n largest contributor is
      # greater than 85% of the extrapolated aggregated value of that cell (Y)
      if(sum(wmaxr)<=nlarge*1.01 && sum(wmax*xmax)>plim*Y){
        # No need to continue the loop if dominance is already TRUE 
        #       print(paste(ifglldat$ehimgid[1], Y, ifglldat$gridvar[1]))
        return(dominance = TRUE)
      } # JON: Not including "else dominance = FALSE" here, it is not necessary. 
      #      Also, there might be a case where dominance = TRUE for nc = 1, but not for nc = 2
    }
  } else {
    
    weight = ifglldat$weight
    dat = ifglldat$gridvar
    wmax = 0
    xmax = 0
    for (ii in 1:length(dat)) {
      lw = min(weight[ii], nlarge-wmax)
      xmax = xmax + dat[ii]*lw
      wmax = wmax + lw
      if (wmax > 0.999*nlarge) break # just to avoid potential numerical issues (FAQ 7.31) from sumWeight >= nlarge
    }
    if(sum(wmax) * 0.999 <=nlarge && xmax > plim*Y$total){
      dominance = TRUE
    } else dominance = FALSE
  }
  dominance
}




