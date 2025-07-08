confid = function(himg, ifgdat, vars, countFeatureOrTotal, mincount, nlarge, plim, domEstat, 
                  checkDominance, checkPpercent,checkReliability, reliabilitySplit, pPercent, userfun, verbose, ...) {
  #  To avoid R CMD check notes
  weight = data = dominance = NULL
  
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
    if (reliabilitySplit & (dim(ifgdat)[1] > 50000 | dim(himg)[1] > 1000)) {
      if (is.logical(reliabilitySplit)) rsplit = dim(ifgdat)[1] %/% 30000
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






mrg_varestim <- function(x, var, strat, PSU, weight, split, verbose, pseudoreg, nhimg){
  t0 = proc.time()[3]
  ID = n = hld = w_sum = wdiff = NULL
  nx = dim(x)[1]
  himgids = unique(x$himgid)
  if (inherits(x, "sf")) x = st_drop_geometry(x)
  if (!missing(PSU)) x$ID = x[[PSU]]  
  icor = 0
  if (missing(strat) || is.null(strat)) {
    x$strat = 1
  } else {
    if (!strat %in% names(x)) stop(paste(strat, "is missing from from the data.frame"))
    x$strat = x[[strat]]  
  }
  
  # Check if any grid cells only have unit value weights
  tt = x %>% group_by(himgid) %>% 
    summarise(wdiff = sum(abs(.data[[weight]] - 1), na.rm = T)) %>% filter(wdiff < 0.1)
  if (dim(tt)[1] > 0) {
    x = x[!(x$himgid %in% tt$himgid), ]
  }
  t1 = proc.time()[3]
  if (verbose > 1) cat("Varestim - finished preprocessing - t= ", round(t1-t0,2), "secs \n")
  if (dim(x)[1] == 0) {
    out_var = data.frame(himgid = himgids, rse = 0)
  } else {
    if (split == 1) {
      df = x
      out_var = vardom(dataset = df, Y= var, H = "strat",
                       PSU = "ID",
                       w_final = weight, Dom = "himgid")$all_result
      t21 = proc.time()[3]
      if (verbose > 1) cat("varestim - finished single split vardom estimation in ", round(t21-t1, 2), "secs \n")
    } else {
      himgid <- unique(x$himgid)    
      #' @importFrom dplyr left_join mutate ungroup group_by distinct case_when n
      #' @importFrom sjmisc split_var 
      df_cl <- left_join(x, data.frame(himgid = himgid, cluster = split_var(himgid, n = split)), by = "himgid")
      out_var <- NULL
      t20 = proc.time()[3]
      if (verbose > 1) cat("varestim - will split reliability calcs in ", split, "cases\n")
      for (isp in 1:split){
        t21 = proc.time()[3]
        df = x[which(df_cl$cluster == isp),]
        if (verbose) print(paste("reliabilitySplit: ", isp, 
                                 "- Number of records: ", paste(dim(df)[1], 
                                                                " - Number of unique IDs: ", length(unique(df$himgid)))))
        
        icor = icor + 1
        t <- df %>% group_by(strat) %>% 
          summarise(hld=n(), w_sum = sum(.data[[weight]], na.rm = T)) %>% 
          filter(hld == 1 & w_sum > 1) %>% ungroup
        if (!is.null(pseudoreg)) pcor = as.numeric(as.factor(df[[pseudoreg]])) else pcor = 0
        if (dim(t)[1] > 0){
          h_st <- t %>% distinct(strat) %>% pull()
          df <- df %>% mutate(strat = case_when(strat %in% c(h_st)~(99999-pcor),
                                                T ~ strat))}
        #' @importFrom vardpoor vardom
        t22 = proc.time()[3]
        if (verbose > 1) cat("varestim - ready to call vardom for split ", isp, "of", split, "time:", round(t22-t21,2), "secs\n")
        est <- vardom(dataset = df, Y= var, H = "strat",
                      PSU = "ID",
                      w_final = weight, Dom = "himgid")
        t23 = proc.time()[3]
        if (verbose > 1) cat("varestim - ready to call vardom for split ", isp, "of", split, "time:", round(t23-t22,2), "secs\n")
        out_var<-rbind(out_var, est$all_result)
      }
    }
    if (verbose) cat("varestim - finished looping in totally ", round(proc.time()[3]-t1,2), "secs \n")
    if (icor > 0) {
      cat(icor, "of the subsets included strata with only one record. \n",
          "You might want to check the strata or consider a lower value for reliabilitySplit.")
    }
    out_var$himgid = as.numeric(out_var$himgid)
    if (dim(tt)[1] > 0) out_var = rbind(out_var[,c("himgid", "rse")], data.frame(himgid = tt$himgid, rse = 0))
    out_var = out_var[order(out_var$himgid),]
  }
  if (!missing(nhimg) && nhimg != dim(out_var)[1]) stop("The dimension of out_var does not match the dimension of himg")
  if (sum(duplicated(out_var$himgid)) > 0) stop("There are duplicated himgids in out_var")
  if (verbose > 1) cat("Finished mrg_varestim")
  out_var$rse
}


pPercentRule = function(ifglldat, pPercent) {
  # Prevent R CMD check notes
  gridvar = weight = NULL
  if (pPercent > 1) pPercent = pPercent/100
  Y <- sum(ifglldat$gridvar*ifglldat$weight)
  ifglldat <- ifglldat[order(ifglldat$gridvar, ifglldat$weight, decreasing = TRUE),] %>%
    mutate(gridvarTotal = gridvar * weight)
  ysec = ifglldat$gridvarTotal[2]
  (Y-ifglldat$gridvarTotal[2]-ifglldat$gridvarTotal[1])/ifglldat$gridvarTotal[1] < pPercent
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
