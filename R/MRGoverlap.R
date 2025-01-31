#' Function that finds and merges overlapping grid cells in a multi-resolution grid
#' The need for this function comes from an error in the gridding process, 
#' and it can be seen as symptom solving rather than solving the issue. 
#' This function will either just show the problematic grid cells or remove the overlaps.
#' 
#' 
#' @eval MRGparam("himg")
#' @eval MRGparam("vars")
#' @eval MRGparam("himg2")
#' @eval MRGparam("action")
#' 
#' @details
#' A multi-resoluion grids should not have overlapping grid cells, by definition.  
#' However, this could happen through stitching different grids together. 
#' Although this should rather have been taken care of during the gridding
#' process, this is not always possible to redo for an end-user. 
#' 
#' This function can first of all be used to identify and show the 
#' overlapping grid cells. It can also be used to create a valid multi-resolution grid for these cases,
#' as long as the grd cells all have the same base grid (i.e. no overlapping grid cells are 
#' partly overlapping, whether the grid cells have the same size or not). 
#' However, there will be some additional errors introduced in this process. 
#' 
#' Except for \code{action = "none"} and \code{action = "replace"}, the function will 
#' try to create a valid grid based on the values in the grid. 
#' This means that it has to find a sensible value for the merged grid cells.
#' Frequently one of the grid cells will have an NA value, meaning that it has been 
#' suppressed. This means that there are observations in the grid cell, but relatively
#' few. If the overlapping grid cell has a value, this is most likely larger, it is
#' non-confidential. The default action is therefore to sum values, but ignore NA-values
#' unless both are NA or one is NA and the other is 0. The last case means that
#' one of the grid cells have a confidential number of records, whereas the 
#' other one has zero records. The total is then a confidential number of records.
#' 
#' if \code{action = "restart"}, there must be a second grid which includes the updated values of 
#' the overlapping grid cells. This could typically happen if the data set is too large
#' to be processed as a single batch. There could then be overlapping grid cells on the border 
#' between different batches. Instead of reprocessing the entire grid, it is possible to
#' reprocess the border regions (in one more more batches, as long as they are not overlapping).
#' The grid cells from the border regions have to be passed as \code{himg2}.
#' 
#' The function uses \code{\link[sf]{st_join}} to check for overlaps. If the grid has
#' a very high number of grid cells (a few tens of thousands), this process can be 
#' rather slow. In that case, it might be better to check parts of the grid separately.
#' 
#' @examples
#' \donttest{
#' library(sf)
#' library(giscoR)
#' 
#' # These are SYNTHETIC agricultural FSS data 
#' data(ifs_dk) # Census data
#' 
#' # Create spatial data
#' ifg = fssgeo(ifs_dk, locAdj = "LL")
#' 
#' nuts2 = gisco_get_nuts(nuts_level = 2)
#' nuts2 = nuts2 %>% filter(CNTR_CODE == "DK") %>% st_transform(crs = st_crs(ifg))
#' 
#' ifg$NUTS2 = st_join(ifg, nuts2, join = st_nearest_feature)$NUTS_ID
#' ress = c(1,5,10,20,40, 80, 160)*1000
#' # Create regular grid of the variables, for three regions
#' ifl = gridData(ifg[ifg$NUTS2 %in% c("DK03", "DK04", "DK05"),], vars = c("UAA"), res = ress)
#' ifl3 = gridData(ifg[ifg$NUTS2 == "DK03",], vars = c("UAA"), res = ress)
#' ifl4 = gridData(ifg[ifg$NUTS2 == "DK04",], vars = c("UAA"), res = ress)
#' ifl5 = gridData(ifg[ifg$NUTS2 == "DK05",], vars = c("UAA"), res = ress)
#'
#' # Create the different multi-resolution grids for different nuts regions
#' himg3 = multiResGrid(ifl3, vars = "UAA", ifg = ifg[ifg$NUTS2 == "DK03",], suppresslim = 0.02)
#' himg4 = multiResGrid(ifl4, vars = "UAA", ifg = ifg[ifg$NUTS2 == "DK04",], suppresslim = 0.02)
#' himg5 = multiResGrid(ifl5, vars = "UAA", ifg = ifg[ifg$NUTS2 == "DK05",], suppresslim = 0.02)
#' 
#' # Bind them together and create new consecutive IDs for the grid cells
#' himg = rbind(himg3, himg4, himg5)
#' himg$ID = 1:dim(himg)[1]
#' 
#' # Find the overlapping grid cells, and show some examples.
#' himgd = MRGoverlap(himg, action = "none")
#' dim(himgd)
#' himgd[himgd$ID.y %in% 932:940,]
#' 
#' # Remove overlapping grid cells
#' himgnew = MRGoverlap(himg, action = "sum")
#' 
#' # Check that there are no more overlaping grid cells
#' himgd2 = MRGoverlap(himgnew, action = "none")
#' himgd2
#' 
#' 
#' # Create a new multi-resolution grid which has the correct grid cells
#' # at the border. In this example, the region of interest is so small that
#' # it is difficult to reprocess just the border grid cells, so 
#' # we make a new complete grid
#' 
#' himg1 =  multiResGrid(ifl, vars = "UAA", ifg = ifg[ifg$NUTS2 %in% c("DK03", "DK04", "DK05"),], suppresslim = 0.02)
#' himgnew2 = MRGoverlap(himg, himg2 = himg1, action = "replace")
#' himgd12 = MRGoverlap(himgnew2, action = "none")
#' himgd12
#' 
#' 
#' }
#' 
#' 
#' @export
MRGoverlap = function(himg, vars, himg2, action = "sum") {
  ID.y = ID.x = NULL # Avoid lacking visible binding
  if (!"ID" %in% names(himg)) himg$ID = 1:dim(himg)[1]  
  if (!"res" %in% names(himg)) himg$res = sqrt(st_area(himg))
  if (missing(vars)) vars = getVars(himg, incCount = TRUE)
  units(himg$res) = NULL
  dups1 = which(duplicated(himg$geometry))
  dups2 = which(duplicated(himg$geometry, fromLast = TRUE))
  hdup = himg[unique(c(dups1, dups2)),]  
  #' @importFrom sf st_join
  hdup = st_join(hdup, hdup, join = st_within)
  hdup = hdup[hdup$ID.x < hdup$ID.y, ]
  
  datjs = NULL
  ress = unique(himg$res)
  for (isep in 2:length(ress)) {
    sep = 0.5*(ress[isep-1] + ress[isep])
    dats = himg[himg$res < sep, ]
    datl = himg[himg$res > sep, ]
    datj = st_join(dats, datl, join = st_within)
    datj = datj[!is.na(datj$ID.y),]
    if (dim(datj)[1] > 0) datjs = rbind(datjs, datj)
  }
  if (dim(hdup)[1] > 0) datjs = rbind(datjs, hdup)  
  datjs = datjs[!duplicated(datjs[, c("ID.x", "ID.y")]),]
  cat("Overlapping grid cells: ", dim(datjs)[1], "\n" )
  if (action == "none") return(datjs)
  if (action == "replace" & !is.null(datjs)) {
    idrem = unique(c(datjs$ID.x, datjs$ID.y))
    hrem = himg[himg$ID %in% idrem, ]
    # Find the rows to merge
    pmerge = st_join(himg2, himg[himg$ID %in% idrem,], join = st_within) %>% 
#      filter(!is.na(ID.x), !is.na(ID.y)) %>%
      filter(!is.na(ID.y)) %>%
      filter(!duplicated(ID.x))
    idadd = unique(pmerge$ID.x)
    himg = himg[!himg$ID %in% idrem,]
    hadd = himg2[himg2$ID %in% idadd, names(himg2) %in% names(himg)]
    himg = rbind(himg[,names(himg) %in% names(himg2)], hadd)
  } else {
  datjs = st_drop_geometry(datjs)
  hrem = NULL
  for (idub in 1:dim(datjs)[1]) {
    id1 = datjs$ID.x[idub]
    id2 = datjs$ID.y[idub]
    hid1 = which(himg$ID == id1)
    hid2 = which(himg$ID == id2)
    for (ivar in 1:length(vars)) {
      v1 = datjs[idub, paste0(vars[ivar], ".x")]
      v2 = datjs[idub, paste0(vars[ivar], ".y")]
      himg[hid2, vars[ivar]] = newval(v1, v2, action)
      wn = paste0("weight_", vars[ivar])
      if (length(grep(wn, names(himg))) > 0) {
        v1 = datjs[idub, paste0(wn, ".x")]
        v2 = datjs[idub, paste0(wn, ".y")]
        himg[hid2, wn] = newval(v1, v2, action)
      }
      if (FALSE) {
      if ("count" %in% names(himg)) {
        v1 = datjs[idub, "count.x"]
        v2 = datjs[idub, "count.y"]
        himg[hid2, "count"] = newval(v1, v2, action)
      }
      if ("countw" %in% names(himg)) {
        v1 = datjs[idub, "countw.x"]
        v2 = datjs[idub, "countw.y"]
        himg[hid2, "countw"] = newval(v1, v2, action)
      }
      }
    }
    hrem = c(hrem, hid1)    
  }
  if (length(hrem) > 0) himg = himg[-hrem, ]
  }
  himg
}

newval = function(v1, v2, action) {
  if (action %in% c("sum", "avg")) {
     if ((is.na(v1) & is.na(v2) | (is.na(v1) & v2 == 0) | (is.na(v2) & v1 == 0))) {
      nv = NA
    } else if (is.na(v1) & v2 > 0) {
      nv = v2
    } else if (is.na(v2) & v1 > 0) {
    nv = v1
    } else if (!is.na(v1) & !is.na(v2)) {
    nv = v1 + v2
    if (action == "avg") nv = nv/2
    } else stop("what happened here 1")
  } else if (action %in% c("sumna", "avgna")) {
    nv = v1 + v2
    if (action == "avgna") nv = nv/2
  } else if (action != "none") stop(cat(action, " is not a possible action parameter \n"))
nv
}
