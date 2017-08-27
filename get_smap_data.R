# library(smapr) # had issues to do with validating username/pass so downlaoded source code to step through
source('~/git/smapr/R/download_smap.R')
source('~/git/smapr/R/extract_smap.R')
source('~/git/smapr/R/find_smap.R')
source('~/git/smapr/R/list_smap.R')
source('~/git/smapr/R/smapr-package.R')
source('~/git/smapr/R/zzz.R')
library("httr")
library("rappdirs")
library("raster")
library("rgdal")
library("rhdf5")
library("xml2")
library("rvest")
library("smapr")

# to download data, need to first set credentials
# see below websites for description of data
# http://nsidc.org/data/docs/daac/smap/sp_l4_sm/data-fields.html#
# https://nsidc.org/data/smap/smap-data.html

if(download.data<-FALSE){
    Sys.setenv(ed_un = 'jamesmaino', ed_pw = 'Jmai0148')
    available_data <- find_smap(id = "SPL3SMAP", date = "2015-05-25", version = 3)
    str(available_data)
    downloads <- download_smap(available_data, directory = "E:\\satellite\\SMAP")
    str(downloads)
    list_smap(downloads, all = FALSE)
    sm_raster <- extract_smap(downloads, "Soil_Moisture_Retrieval_Data/soil_moisture")
    plot(sm_raster, main = "Level 3 soil moisture: May 25, 2015")
}

# takes a while to download so load previously downloaded SPL4SMAP data
# first need to create a data.frame pointing to data that can be read by smapr functions
local_dir<-"E:\\satellite\\SMAP"
dir <- "SPL4SMGP.003/2015.05.25/"
date <- "2015-05-25"
name<-c("SMAP_L4_SM_gph_20150525T013000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T043000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T073000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T103000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T133000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T163000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T193000_Vv3030_001",
        "SMAP_L4_SM_gph_20150525T223000_Vv3030_001")
myfile<-data.frame(name=name, date=date, dir=dir, local_dir=local_dir)
# only need name and local_dir to retrieve files
myfile<-data.frame(name=name, local_dir=local_dir)

# to see contents of SMAP data files use 'list_smap'
list_smap(myfile[1,], all = TRUE)
sm_raster2 <- extract_smap(myfile[2,], "Geophysical_Data/sm_surface", in_memory = TRUE)
sm_raster <- extract_smap(myfile[1,], "cell_lon", in_memory = TRUE)
transformTo <- function(r1){
    ### 0.25*0.25 degree resolution and extent -180, 180, -90, 90
    r=raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
             nrows=1624,ncols=3856,crs="+init=epsg:4326")
    projectRaster(r1,r)
}
rr<-transformTo(sm_raster)
plot(rr, ylim=c(-90,90), xlim=c(-180,180), main = names(rr))
object.size(sm_raster)/10^6
plot(sm_raster)

# crop and plot data for Australia
aus_ex<-extent(c(1.05e7, 15e6, -5.5e6,-1e6)) # australian extent
qld_ex<-extent(c(13315107, 15e6, -3547395,-1e6)) # qld extent
plot(crop(sm_raster, qld_ex))

# plot surface (0-5cm) soil moisture for Melbourne
sm_raster <- extract_smap(myfile, "Geophysical_Data/sm_surface")#, in_memory = TRUE)
loc<-dismo::geocode('Melbourne, VIC',oneRecord = TRUE)
xy<-data.frame(lon=loc$lon, lat=loc$lat)
coordinates(xy) <- c("lon", "lat")
proj4string(xy)<- CRS("+init=epsg:4326") # WGS 84
spTransform(xy, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(1.5 + 3*0:7 -10, raster::extract(sm_raster, xy))


