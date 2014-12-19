## This script calculate richness at different spatial scales from herbarium records.
## It reports raw richness, number of records, and rarefied richness.
## Data will be used to evaluate spatial variation in species area relationships.

############################
## Load packages

library(sp)
library(rgdal)
library(raster)

options(stringsAsFactors=F)

############################
## Functions

####################################
## Access parameters from command when script called

# Get arguments passed from R CMD BATCH
args = commandArgs(trailingOnly = F)
myarg = args[length(args)] # my argument is at the end
grain = as.numeric(sub("-","",myarg))
print(grain)

#grain =  100 # in km

#####################################
## Read in Data

# Sampling grid (produced previously in make_sample_grid.R)
mygrid_ll = readOGR('NAm_sample_grids','NAm_sample_grid_100km')

# Environmental Data (from WorldClim, 10 arc-min resolution)
load('bioclim.Rdata')

# DEM (from WorldClim, 20 arc-sec (1km) resolution)
elev = raster('alt.bil')


#####################################
## Calculate number of records, richness, and rarefied richness for each grid point

# Subset WorldClim data to variables of interest
my.bio = c(1,3,12,15) # Use mean annual temp, isothermality, annual precip, and precip seasonality from WorldClim
bio.stack = stack(bio.stack)
bio.stack = subset(bio.stack, my.bio)

# Convert temperature to kelvin
bio.stack[[1]] = bio.stack[[1]]+273.15

# Calculate mean and variance of climate variables within buffer
clim_mean = extract(bio.stack, mygrid_ll, buffer=grain*1000, fun=mean, na.rm=T)
clim_var = extract(bio.stack, mygrid_ll, buffer=grain*1000, fun=var, na.rm=T)
clim_cv = sqrt(clim_var)/clim_mean
colnames(clim_mean) = paste(colnames(clim_mean), grain, 'mean', sep='_')
colnames(clim_var) = paste(colnames(clim_var), grain, 'var', sep='_')
colnames(clim_cv) = paste(colnames(clim_cv), grain, 'cv', sep='_')

# Calculate mean and variance of elevation within buffer
elev_mean = extract(elev, mygrid_ll, buffer=grain*1000, fun=mean, na.rm=T)
elev_var = extract(elev, mygrid_ll, buffer=grain*1000, fun=var, na.rm=T)
elev_cv = sqrt(elev_var)/elev_mean
elev_mean = data.frame(elev_mean); colnames(elev_mean) = paste('elev', grain, 'mean', sep='_')
elev_var = data.frame(elev_var); colnames(elev_var) = paste('elev', grain, 'var', sep='_')
elev_cv = data.frame(elev_cv); colnames(elev_cv) = paste('elev', grain, 'cv', sep='_')

# Combine data
env_data = data.frame(coordinates(mygrid_ll), clim_mean, elev_mean, clim_var, elev_var, clim_cv, elev_cv)
names(env_data)[1:2] = c('Lon','Lat')

###############################
## Save data
file_name = paste('env_reg_', grain, 'km.csv', sep='')

write.csv(env_data, file_name, row.names=F)









