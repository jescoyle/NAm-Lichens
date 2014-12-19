## This script calculate richness at different spatial scales from herbarium records.
## It reports raw richness, number of records, and rarefied richness.
## Data will be used to evaluate spatial variation in species area relationships.

############################
## Load packages

library(sp)
library(rgdal)
library(vegan) # rarefy
library(maptools) # unionSpatialPolygons
library(stringr) # str_trim

options(stringsAsFactors=F)

############################
## Functions

# A function that finds all records within a given distance buffer from a points
find_recs = function(point, recs, buff){
	dists = spDistsN1(recs, point, longlat=T) # In km
	userecs = recs[dists<buff,]

	return(userecs)
}

# A function that returns a vector counting the number of records for each species binomial
tab_species = function(x){
	if(length(x)>0){

	new_x = get.latinbinom(x)
	splitnames = strsplit(new_x, " ")
	justgenus = new_x[as.numeric(lapply(splitnames, length))==1]
	allspecies = new_x[as.numeric(lapply(splitnames, length))>1]

	if(length(allspecies)>0){
		spgenera = unique(sapply(strsplit(allspecies, " "), function(y) y[1]))
		keepgenera = justgenus[!(justgenus %in% spgenera)]
		allspecies = c(allspecies, keepgenera)
	} else {
		allspecies = justgenus
	}

	as.matrix(table(allspecies))
	} else NULL	
}

# A functions that returns Genus species for a scientificName that may contain subspecific taxonomic info
get.latinbinom = function(x){
	namelist = strsplit(str_trim(x), ' ')
	genus = str_trim(sapply(namelist, function(y) y[1]))
	species = str_trim(sapply(namelist, function(y) ifelse(is.na(y[2]), '',y[2])))
	str_trim(paste(genus, species, sep=' '))
}

####################################
## Access parameters from command when script called

# Get arguments passed from R CMD BATCH
args = commandArgs(trailingOnly = F)
myarg = args[length(args)] # my argument is at the end
grain = as.numeric(sub("-","",myarg))
print(grain)

#grain =  100 # in km
nsamp =  c(20, 50, 100, 200, 500, 1000, 1500) # a vector

#####################################
## Read in Data

# CNALH records in Parmeliaceae and Physciaceae
parmphys = read.csv('Parm_Phys_records_2014-09-20.csv')

# Convert to spatial data
parmphys_sp = parmphys
coordinates(parmphys_sp) = c('decimalLongitude','decimalLatitude')
proj4string(parmphys_sp) = CRS("+proj=longlat")

# Sampling grid (produced previously in make_sample_grid.R)
mygrid_ll = readOGR('NAm_sample_grids','NAm_sample_grid_100km')

#####################################
## Calculate number of records, richness, and rarefied richness for each grid point

# Define data frame where going to store data
mygrid_data = mygrid_ll

# Calculate number of records and raw richness
nobs = sapply(1:length(mygrid_ll), function(i){
	userecs = find_recs(mygrid_ll[i,], parmphys_sp, grain)
	
	if(nrow(userecs)>0){
		parm = userecs[userecs$family=='Parmeliaceae',]
		phys = userecs[userecs$family=='Physciaceae',]

		N_parm = nrow(parm)
		N_phys = nrow(phys)

		S_parm = length(tab_species(parm$scientificName))
		S_phys = length(tab_species(phys$scientificName))
	} else {
		N_parm = 0
		N_phys = 0
		S_parm = 0
		S_phys = 0
	}

	c(N_parm, N_phys, S_parm, S_phys)
})
nobs_df = data.frame(t(nobs))
names(nobs_df) = c('N_parm','N_phys','S_parm','S_phys')

mygrid_data = cbind(mygrid_data, nobs_df)
names(mygrid_data)[2:3] = c('Lon','Lat')

# Calculate rarefied richness
rich_rare = data.frame()

for(i in 1:length(mygrid_ll)){
	userecs = find_recs(mygrid_ll[i,], parmphys_sp, grain) # 500km radius

	if(nrow(userecs)>0){

	parm = userecs[userecs$family=='Parmeliaceae',]
	phys = userecs[userecs$family=='Physciaceae',]

	sapply(nsamp, function(n){
		if(nrow(parm)>=n){
			comm = tab_species(parm$scientificName)
			as.numeric(rarefy(comm, n))
		} else {
			NA
		}
	}) ->  rich_parm

	sapply(nsamp, function(n){
		if(nrow(phys)>=n){
			comm = tab_species(phys$scientificName)
			as.numeric(rarefy(comm, n))
		} else {
			rich = NA
		}
	}) ->  rich_phys

	} else {
		rich_parm = rep(0, length(nsamp))
		rich_phys = rep(0, length(nsamp))
	}
	
	rich_rare = rbind(rich_rare, c(rich_parm, rich_phys))

}

names(rich_rare) = c(paste('S',nsamp,'_parm', sep=''),paste('S',nsamp,'_phys', sep=''))

# Add to mygrid_data
mygrid_data = cbind(mygrid_data, rich_rare)


###############################
## Save data
file_name = paste('parmphys_richgrid_', grain, 'km.csv', sep='')

write.csv(mygrid_data, file_name, row.names=F)









