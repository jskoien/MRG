
MRGparam = function(par) {

  ret = switch(par,
               MRGobject = "@param MRGobject An object including all
                 the necessary variables and parameters, from a call to createMRGobject",
               gdl = "@param gdl A list of gridded data - with different resolutions",
              himg = "@param himg The grid resulting from a call to multiResGrid",
               ifs = "@param ifs A data.frame or tibble with the locations and the data of the survey or census data" ,
               crsOut = "@param crsOut The coordinate reference system (crs) to be used ",
               ifg = "@param ifg Either a data.frame or tibble or sf-object with the locations and the data of the survey or census data,
                          or a list of such objects." ,
               ress = "@param ress A vector with the different resolutions", 
              res = "@param res A resolution or a vector with the different resolutions", 
              lnames = "@param lnames Names for the different surveys or censuses if ifg is a list. 
                         Typically it could be survey years",
               geovar = "@param geovar Name of geodata variable in the objects. Must me the same for all of the surveys/censuses, if 
                              the data sets are not submitted as sf-objects",
               vars = "@param vars Variable(s) of interest that should be aggregated (necessary when ifg is
                        used for individual farm specific anonymization rules)",
               weights = "@param weights Extrapolation factor(s) (weights) wi of unit i in the sample of units nc 
                      falling into
                               a specific cell c. Weights are used for disclosure control measures. 
                               A weight of 1 will be used if missing.
                             If only one weight is given, it will be used for all variables. If the length is more than one,
                             the length has to be equal to the number of variables. If the same weight is used for several variables,
                             it must be repeated in the weights-vector",
               mincount = "@param mincount The minimum number of farms for a grid cell (threshold rule)",
              countFeatureOrTotal = "@param countFeatureOrTotal Should the frequency limit be applied on records with a positive
                               value for a certain feature, or on all records, independent of value of feature", 
              minpos = "@param minpos Minimum number of positive values for a variable", 
               nlarge = "@param nlarge Parameter to be used if the nlarge(st) farms should count for maximum plim percent of
                        the total value for the variable in the grid cell (see details of \\code{\\link{gridData}})",
               plim = "@param plim See nlarge",
               verbose = "@param verbose indicates if some extra output should be printed",
               nclus = "@param nclus Number of clusters to use for parallel processing. No parallelization is used
                           for \\code{nclus = 1}.",
               clusType = "@param clusType The type of cluster; see \\code{\\link[parallel]{makeCluster}} for more details.
                        The default of makeCluster is used if type is missing or NA",
               domEstat = "@param domEstat Should the dominance rule be applied as in the IFS handbook (TRUE), where 
                         the weights are rounded before finding the  first nlarge contributors, or should 
                        it be the first nlarge contributors*weight, where also fractions are considered (FALSE)?", 
               consistencyCheck = "@param consistencyCheck logical; whether consistency between the gridded values and
                         the similar values from ifg should be checked. The gridded value is derived
                         from rasterize and the second one from st_join. The two methods can in some 
                         cases treat border cases between grid cells differently.",       
               outfile = "@param  outfile File to direct the output in case of parallel processing, 
                         see \\code{\\link[parallel]{makeCluster}} for more details.",
               splitlim = "@param splitlim For large dataset - split the data set in batches of more or less splitlim size",
               checkDominance = "@param checkDominance Logical - should the dominance rule be applied?",
               checkReliability = "@param checkReliability Logical - should the prediction variance be checked, and used for the aggregation?
                         This considerably increases computation time",
               userfun = "@param userfun This gives the possibility to add a user defined function with additional confidentiality rules which 
                             the grid cell has to pass", 
              fargs = "@param fargs The name of the necessary variables of userfun", 
              strat = "@param strat Column name defining the strata for stratified sampling, used if checkReliability is TRUE",
               confrules = "@param confrules Should the frequency rule (number of holdings) refer to the number of holdings with 
                           a value of the individual vars above zero (\"individual\") or the total number of holdings in 
                           the data set (\"total\")? ",
               suppresslim = "@param suppresslim Parameter that can be used to avoid that almost empty grid cells are merged with cells 
                            with considerably higher number of observations. The value is a minimum percentage of the total
                            potential new cell for a grid cell to be aggregated.",  
               sumsmall = "@param sumsmall Logical; should the suppresslimSum value be applied on the sum of
                            small grid cells within the lower resolution grid cell?
                            Note that different combinations of suppreslim and suppreslimSum values 
                            might not give completely intuitive results.For instance, if both are equal, then
                            a higher value can lead to more grid cells being left unaggregated for smaller grid sizes, leading
                            to aggregation for a large grid cell",
               suppresslimSum = "@param suppresslimSum Parameter similar to suppreslim, but affecting the total
                              of grid cells to be suppressed",
               plotIntermediate = "@param plotIntermediate Logical or number - make a simple plot showing which grid cells have already 
                           passed the frequency rule. plotintermediate = TRUE, the function will wait 5 seconds after plotting 
                           before continuing, otherwise it will wait plotintermediate seconds.",
               addIntermediate = "@param addIntermediate Logical; will add a list of all intermediate himgs
                           and lohs (overlay of himg and the lower resolution grid) as an attribute to
                           the object to be returned",
               reliabilitySplit = "@param reliabilitySplit Logical or number - parameter to be used in
                            calculation of the reliability (if checkReliability = TRUE), a limit for when
                            the reliability calculation should be done in batches. 
                              If TRUE it will use the default value (50,000), otherwise a user value. If FALSE, the data set will not
                               be split independent on the size.", 
               locAdj = "@param locAdj parameter to adjust the coordinates if they are exactly on the borders between grid cells. The values
                          can either be FALSE, or \"jitter\" (adding a small random value to the coordinates, essentially spreading
                          them randomly around the real location), \"UR\", \"UL\", \"LR\" or \"LL\", to describe which corner of the grid 
                          cell the location belong (upper right, upper left, lower right or lower left).",
              remZeroes = "@param remZeroes Set to TRUE if the gridding should only be done on a 
                            subset of the data set, for which the value(s) of the variable(s) of interest
                            is larger than zero. If there is more than one variable of interest, the subset
                            will include all records where at least one of the variables are above zero.",
              postProcess = "@param postProcess Logical; should the postprocessing be done as part
                             of creation of the multiresolution grid (TRUE), or be done in a separate 
                             step afterwards (FALSE). The second option is useful when wanting
                             to check the confidential grid cells of the final map",
              remCols = "@param remCols Logical; Should intermediate columns be removed? Can be set
                          to FALSE for further analyses",
              rounding = "@param rounding either logical (FALSE) or an integer indicating the number 
                      of decimal places 
                       to be used. Negative values are allowed (such as the default
                       value rounding to the closest 10). See also the details
                       for \\code{digits} in \\code{\\link{round}}.",
              ellipsis = "@param ... Possible arguments to userfun or other internal functions"
              
)
ret
}
               
