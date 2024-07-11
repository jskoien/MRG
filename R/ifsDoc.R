#' Test data sets for the multiResGrid package
#' 
#' 
#' Synthetic data set of Danish farming data, similar to the 
#' structure of the real Farm Structure Survey (FSS) data. It contains more than 
#' 37,000 synthetic records - generated in a way that should
#' replicate the structure and the distribution of real data,
#' but where the individual data are different from the real data.
#' 
#' The variables are as follows:
#' @format A data frame with 37,088 rows and 14 variables
#' \itemize{
#'   \item COUNTRY The name of the country
#'   \item YEAR The year of the survey data
#'   \item ID_SYNTH Unique ID of the record
#'   \item FARMTYPE Farm typology. Farms are classified into different types according to their dominant activity and standard output value (proxy for farm income). For further information see https://ec.europa.eu/eurostat/statistics-explained/index.php?title=Glossary:Farm_typology
#'   \item HLD_FEF Not used. Farm is included in frame extension (HLD_FEF=1) or main frame (HLD_FEF=0)
#'   \item REGIONS NUTS2 region
#'   \item GEO_LCT The geolocation in typical FSS-format, including both country, CRS and xy coordinates
#'   \item EXT_CORE The extrapolation weights for core data (1 in this data set)
#'   \item STRA_ID_CORE Which stratum the record belongs to - only used for the reliability checking
#'   \item UAA The utilized agricultural area of the farm
#'   \item UAAXK0000_ORG The organic utilized agricultural area, excluding kitchen gardens of the farm. UAAXK0000_ORG includes fully certified area and area under conversion 
#'   \item Sample Whether the record should be included as a weighted subsample
#'   \item EXT_MODULE The extrapolation weights for the sample data 
#' }
#' @docType data
#' @keywords datasets
#' @name ifs_dk
#' @usage data(ifs_dk)
#' @format A data frame with 37088 rows and 14 variables
NULL
