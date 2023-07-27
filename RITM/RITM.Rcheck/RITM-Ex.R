pkgname <- "RITM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "RITM-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('RITM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("areaT")
### * areaT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: areaT
### Title: Area Mode Irregular Terrain Modeling
### Aliases: areaT

### ** Examples


ModVar = 3;
deltaH = 90;
tht_m = 100;
rht_m = 10;
dist_km = 20;
TSiteCriteria = 0;
RSiteCriteria = 0;
eps_dielect = 15;
sgm_conductivity = 0.005;
eno_ns_surfref = 301;
frq_mhz = 145;
radio_climate = 1;
pol = 1;      #1 = vert
pctTime = 0.5;
pctLoc = 0.5;
pctConf = 0.9;

areaT(ModVar, deltaH, tht_m, rht_m, dist_km, TSiteCriteria, RSiteCriteria,
eps_dielect, sgm_conductivity, eno_ns_surfref,frq_mhz, radio_climate, pol, pctTime, pctLoc,
pctConf)$dbloss




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("areaT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("point_to_point")
### * point_to_point

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: point_to_point
### Title: Point to Point ITM
### Aliases: point_to_point

### ** Examples


#commented below is an example of how to get an elevation profile in R.
#library(elevatr)
#library(sp)
#set.seed(65.7)
#examp_df <- data.frame(x = runif(3, min = -73, max = -72.5), y = runif(3, min = 42,  max = 43))
#prj_dd <- "+init=EPSG:4326"
#cats <- data.frame(category = c("H", "M", "L"))
#examp_df2 <- data.frame(examp_df, cats)
#examp_sp <- SpatialPoints(examp_df, proj4string = CRS(prj_dd))
#examp_spdf <- SpatialPointsDataFrame(examp_sp, data = cats)
#df_elev_epqs <- get_elev_point(examp_df, prj = prj_dd, src = "epqs")
#Elevation<-df_elev_epqs$elevation

#These are the values returned above:
Elevation<-c(207.81, 198.95, 306.15)

#Build the input list
struct_Input<-list()
struct_Input$Frequency<-120*1000000 #Frequency to calculate loss at (Hz)
struct_Input$Elevation<-Elevation #terrain elevation profile, (list of points) (m)
struct_Input$Resolution<-40000 #terrain input resolution (distance b/t points) (m)
struct_Input$TX_Height<-3 #Transmit antenna height above ground (m)
struct_Input$RX_Height<-100 #Recieve antenna height above ground (m)
struct_Input$eps<-15 #Soil dielectric
struct_Input$sgm<-.005 #Surface conductivity
struct_Input$surfref<-301 #Surface refractivity
struct_Input$Climate<-5
struct_Input$Polarization<-1 #1 is vertical, 0 is horizontal
struct_Input$Confidence<-.95 #confidence for statistical analysis
struct_Input$Reliability<-.95 #Reliability to calculate statistics for (.01 to .99)

point_to_point(struct_Input)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("point_to_point", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
