## This script produces a grid of points across terrestrial North America

## Load packages

library(sp)
library(rgdal)
library(maptools) # unionSpatialPolygons

options(stringsAsFactors=F)

# Map of North America
nam_outline = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
nam_outline = nam_outline[ nam_outline$COUNTRY %in% c('USA','MEXICO','CANADA'), ]

## Create sampling grid

# Convert to equal-area projection
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')
nam_laea = spTransform(nam_outline,CRS(plot_prj))

# Find bounding box and make grid
grid_scale = 100 # in km approximately, since laea preserves area not distance
mybox = bbox(nam_laea)
x = seq(mybox[1,1]+50, mybox[1,2], grid_scale)
y = seq(mybox[2,1]+50, mybox[2,2], grid_scale)
mygrid_laea = expand.grid(x,y)

# Exclude grid points not on land
nam_poly = SpatialPolygons(nam_laea@polygons, proj4string=CRS(plot_prj))
nam_poly1 = unionSpatialPolygons(nam_poly, IDs = rep(1,length(nam_poly)))
mygrid_sp = SpatialPoints(mygrid_laea, proj4string=CRS(plot_prj))
usegrid = over(mygrid_sp, nam_poly1)
mygrid_sp = mygrid_sp[!is.na(usegrid),]

# View grid
#plot(nam_poly1)
#plot(mygrid_sp, add=T)

# Reproject these coordinates to lat/lon
mygrid_ll = spTransform(mygrid_sp, CRS("+proj=longlat +datum=WGS84"))
plot(mygrid_ll)

## Save grid
layer_name = paste('NAm_sample_grid_', grid_scale, 'km', sep='')
writeOGR(mygrid_ll, 'NAm_sample_grids', layer_name, 'ESRI Shapefile')
