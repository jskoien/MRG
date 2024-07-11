library(MRG)
library(sf)
library(giscoR)
#'
# These are SYNTHETIC agricultural FSS data 
data(ifs_dk) # Census data
ifs_weight = ifs_dk %>% dplyr::filter(Sample == 1) # Extract weighted subsample

# Create spatial data
ifg = fssgeo(ifs_dk, locAdj = "LL")
fsg = fssgeo(ifs_weight, locAdj = "LL")
# Read country borders, only used for plotting
borders = gisco_get_nuts(nuts_level = 0)
dkb = borders[borders$CNTR_CODE == "DK",] %>% st_transform(crs = 3035)
#'
# Set the base resolutions, and create a hierarchical list with gridded data
ress = c(1,5,10,20,40, 80, 160)*1000
# Gridding Utilized agricultural area (UAA)
ifl = gridData(ifg, "UAA",res = ress)

# Gridding UAA and organic UAA together
ifl3 = gridData(ifg, vars = c("UAA", "UAAXK0000_ORG"), res = ress)

# Gridding the UAA from the survey - the survey weights are in the column EXT_MODULE
fsl = gridData(fsg,  vars = c("UAA"), weights = "EXT_MODULE",  res = ress)

# Create a multi-resolution grid only with farm number as confidentiality rule, then plot results
himg0 = multiResGrid(ifl, checkReliability = FALSE, suppresslim = 0)

# Create a multi-resolution grid of UAA, also based on the dominance rule (default)
himg1 = multiResGrid(ifl, vars = "UAA", ifg = ifg)


# Create joint multi-resolution grid of organic UAA and total UAA
himg3 = multiResGrid(ifl3, vars = c("UAA", "UAAXK0000_ORG"), ifg = ifg, 
                  checkReliability = FALSE, suppresslim = 0)


# Create joint multi-resolution grid of organic UAA and total UAA
himg4 = multiResGrid(ifl3, vars = c("UAA", "UAAXK0000_ORG"), ifg = ifg, 
                     checkReliability = FALSE, suppresslim = 0.1)

# Create multi-resolution grid of UAA and organic UAA, based on survey data,
# also applying reliability check
# Slow!
himg5 = multiResGrid(fsl, vars = c("UAA"), weights = "EXT_MODULE", ifg = fsg, 
                      strat = "STRA_ID_CORE", checkReliability = TRUE)

summary(himg0)
summary(himg1)
summary(himg3)
summary(himg4)
summary(himg5)

MRGobject = createMRGobject(ifg = ifg, ress = ress, var = "UAA")
himg1 = multiResGrid(MRGobject)
# Parameters can be updated in the object or in the call to multiResGrid
MRGobject$suppresslim = 0.02
himg2 = multiResGrid(MRGobject)
himg3 = multiResGrid(MRGobject, suppresslim = 0.05)
summary(himg1)
summary(himg2)
summary(himg3)


