## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(sf)
library(viridis)
library(ggplot2)
library(giscoR)
library(MRG)

data(ifs_dk)
# Create spatial data
ifg = fssgeo(ifs_dk)

# Read country borders of Denmark, only used for plotting
borders = gisco_get_nuts(nuts_level = 0)
dkb = borders[borders$CNTR_CODE == "DK",] %>% st_transform(crs = 3035)


# Set the base resolutions, and create a hierarchical list with gridded data
# No variable name is needed if we're only looking at the number of holdings
# but we include the utilizied agricultural area, as it will also be used further below.
# The jitter is necessary to avoid border problems, creating grid borders exactly 
# on the 
# locations of the holdings (rounded to closest 1000m)
ress = 1000*2^(1:7)
ifl = list()
ifg = st_jitter(ifg, amount = 50)
ifl = gridData(ifg, vars = "UAA",  res = ress)

# Run the adaptive grid function only with farm number as confidentiality rule, 
# then plot results
himg1 = multiResGrid(ifl, checkValidity = FALSE)


## -----------------------------------------------------------------------------

ggplot() + geom_sf(data = himg1, aes(fill = count, color = count)) +
   scale_fill_viridis( name = "number of farms", trans = "log10") +
   scale_color_viridis( name = "number of farms", trans = "log10") +
   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
   ggtitle("Number of farms for variable pixel size, only number of farms confedentiality") +
   theme_bw()

## -----------------------------------------------------------------------------
 himg2 = multiResGrid(ifl, ifg = ifg, vars = "UAA", checkValidity = FALSE)
ggplot() + geom_sf(data = himg2, aes(fill = count, color = count)) +
   scale_fill_viridis( name = "number of farms", trans = "log10") +
   scale_color_viridis( name = "number of farms", trans = "log10") +
   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
   ggtitle("Number of farms, based on farm numbers and size") +
   theme_bw()

ggplot() + geom_sf(data = himg2, aes(fill = UAA, color = UAA)) +
   scale_fill_viridis( name = "UAA", trans = "log10") +
   scale_color_viridis( name = "UAA", trans = "log10") +
   geom_sf(data = dkb, fill = NA, colour='black', lwd = 1) +
   coord_sf(crs = 3035) +#, xlim = c(2377294, 6400000), ylim = c(1313597, 5628510)) +
   ggtitle("UAA, based on farm numbers and size") +
   theme_bw()


