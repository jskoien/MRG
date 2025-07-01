confid = function(himg, ifgdat, vars, countFeatureOrTotal, mincount, nlarge, plim, domEstat, 
                  checkDominance, checkPpercent,checkReliability, reliabilitySplit, pPercent, userfun, verbose, ...) {
  #  To avoid R CMD check notes
  weight = data = ifg = dominance = NULL
  
  himgdat = st_drop_geometry(himg)
  if (tolower(countFeatureOrTotal) == "feature" & !missing(vars)) {
    ww = himgdat[,names(himgdat) %in% paste0("weight_", vars), drop = FALSE]
    ww = apply(ww, MARGIN = 1, FUN = function(x) if (sum(x > 0))  min(x[x>0]) else 0)
  } else ww = himgdat[, "countw"]
  wf = which(ww > 0 & ww < mincount)
  if (length(wf) > 0)  himg$freq[which(himg$ID %in% himgdat$ID[wf])] = TRUE
  len = ifelse(missing(vars), 1, length(vars))
  ehimgid = ufres = NULL
  for (ivar in 1:len) {
    ifgdatl = NULL
    if ((checkDominance | checkPpercent) & !missing(vars)) {
      ifgdatl <- ifgdat[,c("himgid", paste0("gridvar", ivar), paste0("weight",ivar))] 
      names(ifgdatl) = c("himgid", "gridvar", "weight")
      ifgdatl$ehimgid = ifgdatl$himgid
      ifgdatl = ifgdatl[order(ifgdatl$ehimgid),]
    }
    #' @importFrom tidyr unnest nest 
    #' @importFrom purrr map 
    if (checkDominance & !missing(vars)) {
      if (verbose) cat("Checking dominance \n")
      dom = ifgdatl %>% filter(weight != 0)  %>%
        group_by(ehimgid) %>%  nest() %>%
        mutate(dominance = map(data, ~dominanceRule(., nlarge = nlarge, plim = plim, 
                                                    domEstat = domEstat))) %>%
        unnest(dominance) %>% ungroup %>% select(dominance) %>% pull 
      
      domid = which(dom)
      if (length(domid) > 0) himg$dom[domid] = TRUE
    } 
    if (checkPpercent & !missing(vars)) {
      if (verbose) cat("Checking p-percent rule \n")
      pPercentC = ifgdatl %>% filter(weight != 0)  %>%
        group_by(ehimgid) %>%  nest() %>%
        mutate(pPercentC = map(data, ~pPercentRule(., pPercent = pPercent))) %>%
        unnest(pPercentC) %>% ungroup %>% select(pPercentC) %>% pull 
      
      pPerid = which(pPercentC)
      if (length(pPerid) > 0) himg$pPerc[pPerid] = TRUE
      
    }
    if (!missing(userfun) && is.function(userfun)) {
      if (verbose) cat("Checking userfun \n")
      if (is.null(ifgdatl)) ifgdatl <- ifgdat[,c("himgid", paste0("gridvar", ivar), 
                                                 paste0("weight",ivar))] 
      
      ifgdatl$ehimgid = ifgdatl$himgid
      ifgdatl = ifgdatl[order(ifgdatl$ehimgid),]
      dots = list(...)
      # The next 14 lines are based on an answer on StackOverflow:
      # https://stackoverflow.com/questions/78647845/using-purrrmap-with-a-user-defined-function-how-to-pass-arguments
      # by the user Nir Graham (userid: 11726436)
      # licensed by StackOverlow under CC BY-SA 4.0
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
      if (length(ufid) > 0) himg$ufun[ufid] = TRUE
    }
  }
  
  if (checkReliability) {
    if (verbose) cat("Checking reliability \n")
    rsplit = reliabilitySplit
    if (reliabilitySplit & (dim(ifg)[1] > 50000 | dim(himg)[1] > 1000)) {
      if (is.logical(reliabilitySplit)) rsplit = dim(ifg)[1] %/% 30000
      if (dim(himg)[1] > rsplit*1000) rsplit = dim(himg)[1] %/% 1000
    }
    if (!missing(vars) && !is.null(vars)){
      for (ivar in 1:length(vars)){
        nhimg = dim(himg)[1]
        vestres = mrg_varestim(ifgdat, var = paste0("gridvar", ivar), strat = "strat", PSU = "ID", 
                               weight = paste0("weight", ivar), split = rsplit, pseudoreg = "pseudoreg", 
                               verbose = verbose, nhimg = nhimg)
        himg[,paste0("vres",ivar)] = vestres
      }
    }
    nonvalids = suppressWarnings(which(apply(st_drop_geometry(himg[, grep("vres", names(himg))]), 1, max, na.rm = TRUE) > 0.35))
    himg$reliability[nonvalids] = TRUE
  }
  
  himg$confidential = rowSums(st_drop_geometry(himg[,c("freq", "dom", "pPerc", "ufun", "reliability")])) > 0
  himg
  
}