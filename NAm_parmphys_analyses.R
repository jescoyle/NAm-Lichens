## This script conducts analyses of Parmeliaceae and Physciaceae diversity in North America
## The base data is gridded (rarefied) richness maps at various scales computed by 'calc_regS_across_scales.R' on 11/20/2014

library(sp)
library(rgdal)
library(raster)


############################################################
### Read in Data

working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CNALH Diversity'
data_dir = './Gridded Richness'

setwd(working_dir)

# Read in spatial grid
mygrid_ll = readOGR('NAm_sample_grids','NAm_sample_grid_100km')

# Read in map data as a long data frame
scales  = c(10,22,50,112,250,559,1250)

richness = data.frame()
for(s in scales){
	filename = paste('parmphys_richgrid_',s,'km.csv', sep='')
	this_data = read.csv(paste(data_dir,filename, sep='/'))
	this_data$radius = s
	richness = rbind(richness, this_data)	
}


# Convert df to spatial
richness_sp = richness
coordinates(richness_sp) = c('Lon','Lat'); proj4string(richness_sp) = CRS("+proj=longlat")
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')
richness_sp = spTransform(richness_sp, CRS(plot_prj))

# Read in environmental data
env = data.frame(ID = 1:nrow(mygrid_ll))
for(s in scales[2:5]){
	filename = paste('env_reg_',s,'km.csv', sep='')
	this_data = read.csv(paste(data_dir,filename, sep='/'))
	env = cbind(env, this_data)
}

# North America Outline
nam_outline = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
nam_outline = nam_outline[ nam_outline$COUNTRY %in% c('USA','MEXICO','CANADA'), ]
nam_laea = spTransform(nam_outline,CRS(plot_prj))

#############################################################
### Species Richness Maps

## Decide rarefaction subsample size
sizes = c(20,50,100,200,500,1000,1500)

# count number of pixels with N or fewer species for each rarefaction subsample size
combos = expand.grid(scales, sizes)
counts = array(NA, dim=c(2, length(scales), length(sizes)), 
	dimnames=list(c('parm','phys'), radius=scales, nsamp=sizes))
for(i in 1:nrow(combos)){
	use_data = subset(richness, radius==combos[i,1])
	props = c(sum(use_data$N_parm>=combos[i,2]), sum(use_data$N_phys>=combos[i,2]))/nrow(use_data)
	
	counts[,as.character(combos[i,1]), as.character(combos[i,2])] = props

}


counts['parm',,] # Let's look at maps for nsamp = 100 (all scales) and radius = 250 (all nsamp)
counts['phys',,]

## Make maps

# Color Ramp
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]

# North America Outline
nam_outline = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
nam_outline = nam_outline[ nam_outline$COUNTRY %in% c('USA','MEXICO','CANADA'), ]
nam_laea = spTransform(nam_outline,CRS(plot_prj))

scale_bar = list("SpatialPolygonsRescale", layout.scale.bar(), 
	offset = c(0,0), scale = 1000, fill=c("transparent","black"))


# Map for parmeliaceae, nsamp=100
colcuts = seq(1,73,8)
ncuts=length(colcuts) - 1
mycolramp = colorRampPalette(mycol)(ncuts)

pdf('./Maps/Parmeliaceae richness 100 records.pdf', height=6, width=6)
for(s in scales){
	this_rich = subset(richness_sp, radius==s)

	#trellis.par.set(axis.line=list(col=NA))
	print(

	spplot(this_rich, 'S100_parm',  main='', panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
		sp.text(c(-4500,4500), paste('Radius =',s,'km'))
		
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-'))))
	)

	)
}
dev.off()

# Map for physciaceae, nsamp=100
colcuts = seq(1,73,8)
ncuts=length(colcuts) - 1
mycolramp = colorRampPalette(mycol)(ncuts)

pdf('./Maps/Physciaceae richness 100 records.pdf', height=6, width=6)
for(s in scales){
	this_rich = subset(richness_sp, radius==s)

	#trellis.par.set(axis.line=list(col=NA))
	print(

	spplot(this_rich, 'S100_phys',  main='', panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
		sp.text(c(-4500,4500), paste('Radius =',s,'km'))
		
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-'))))
	)

	)
}
dev.off()

# Map for parmeliaceae, radius = 250km
colcuts = seq(0,210,15)
colcuts[1] = 1
ncuts=length(colcuts) - 1
mycolramp = colorRampPalette(mycol)(ncuts)

this_rich = subset(richness_sp, radius==250)

pdf('./Maps/Parmeliaceae richness 250km radius.pdf', height=6, width=6)
for(n in sizes){
	
	plot_col = paste('S',n,'_parm', sep='')
	print(

	spplot(this_rich, plot_col,  main='', panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
		sp.text(c(-4500,4500), paste('# Records =',n))
		
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-'))))
	)

	)
}
dev.off()

pdf('./Maps/Physciaceae richness 250km radius.pdf', height=6, width=6)
for(n in sizes){
	
	plot_col = paste('S',n,'_phys', sep='')
	print(

	spplot(this_rich, plot_col,  main='', panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
		sp.text(c(-4500,4500), paste('# Records =',n))
		
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-'))))
	)

	)
}
dev.off()


#############################################################
### Species Area Relationships

## Find pixels with non-NA richness at the most scales and nsamps
n = 200

count_data = sapply(unique(richness$ID), function(x){
	this_data = subset(richness, ID==x)

	nscales_parm = sum(this_data[,paste('S',n,'_parm', sep='')]>0, na.rm=T)
	nscales_phys = sum(this_data[,paste('S',n,'_phys', sep='')]>0, na.rm=T)

	c(nscales_parm, nscales_phys)
})

count_data = data.frame(t(count_data))
names(count_data) = c('parm','phys')

hist(count_data$parm)
hist(count_data$phys)

count_data$sum = rowSums(count_data)
count_data[order(count_data$sum),]

#potentials = subset(count_data, sum >= 5)
#pot_sp = subset(richness_sp, ID %in% rownames(potentials))

count_data_sp = count_data
this_rich = subset(richness_sp, radius==250)
coordinates(count_data_sp) = coordinates(this_rich[order(this_rich$ID),])
proj4string(count_data_sp) = CRS(proj4string(this_rich))

colcuts = -1:13
ncuts=length(colcuts)-1
mycolramp = colorRampPalette(mycol)(ncuts)
spplot(count_data_sp[,], 'sum',  main='', panel=function(x,y,subscripts,...){
	panel.pointsplot(x,y,...)
	sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='# scales',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(as.character(0:13)))
)

use_col = mycolramp[cut(count_data_sp$sum, colcuts)]
plot(count_data_sp, col=use_col, pch=15, cex=.7)
identify(coordinates(count_data_sp))

## For five locations: Pacific NW, SW, Maine, Wisconsin, Southern Mexico
# based on choosing pixels from previous map
myregions = data.frame(
	region=c('Appalachians','Sonoran Desert','Pacific Northwest','Upper Midwest','Northeast','Mexico','Alaska','Colorado'),
	ID=c(429,265,1044,940,757,41,2023,611)
)
myregions = myregions[order(myregions$ID),]
Appalachians: 429
Sonoran: 265
PNW: 1044
Midwest: 940
Northeast: 757
Mexico: 41
Alaska: 2023
Colorado: 611


## Plot species area curves
richness$area = pi*(richness$radius^2)

pdf('SAR 8 regions 200 records.pdf', height=5, width=9)
par(mfrow=c(1,2))
plot(0,0, xlim=c(2,7), ylim=c(0,2.5), type='n', las=1, xlab='Log Area',ylab='Log Richness', main='Parmeliaceae')
for(i in 1:8){
	points(log10(S200_parm)~log10(area), data=subset(richness, ID==myregions[i,'ID']), xlim=c(0,7), ylim=c(0,2.5), type='o', col=i, pch=rep(15:18, 2)[i])
}
plot(0,0, xlim=c(2,7), ylim=c(0,2.5), type='n', las=1, xlab='Log Area',ylab='Log Richness', main='Physciaceae')
for(i in 1:5){
	points(log10(S200_phys)~log10(area), data=subset(richness, ID==myregions[i,'ID']), xlim=c(0,7), ylim=c(0,2.5), type='o', col=i, pch=rep(15:18, 2)[i])
}
legend('bottomright', as.character(myregions$region), col=1:8, pch=15:18, lwd=2, lty=1, bty='n')

dev.off()




#############################################################
### Richness Environment Relationships

# Read in previously  made data frames including test/fit dataset partition
data_ll = readOGR('Spatial Model Data','model_data_250km_100_records')
data_sp = spTransform(data_ll, CRS(plot_prj))

## Make dataframe of data for models if not read in above
# Decide on focal scale of analysis and number of records subsampled
focal_scale = 250
focal_n = 100
use_rich = subset(richness, radius==focal_scale)[,c('ID','Lon','Lat',paste('S',focal_n,'_parm',sep=''),paste('S',focal_n,'_phys',sep=''))]
use_env = data.frame(ID=env[,'ID'], env[,grep(focal_scale, colnames(env))])
colnames(use_rich)[4:5] = c('S_parm','S_phys')
colnames(use_env)[2:ncol(use_env)] = sapply(colnames(use_env)[2:ncol(use_env)], function(x){
	this_split = unlist(strsplit(x, '_'))
	paste(this_split[1],this_split[3],sep='_')
})

use_data = merge(use_rich, use_env, by='ID')

# Convert to spatial points data frame
data_ll = SpatialPointsDataFrame(mygrid_ll, use_data)
data_sp = spTransform(data_ll, CRS(plot_prj))

## Partition grid into fit and test data sets from grid where both Parm and Phys have richness

# Determine which sampled points have observed richness for both Parm and Phys
use_obs = apply(use_data[,c('S_parm','S_phys')], 1, function(x) sum(is.na(x))==0)
sum(use_obs) # 954
data_sp$model_obs = use_obs
data_ll$model_obs = use_obs

# Randomly assign points to fit and test data sets
N = 477 # number of points to be assigned to the test data set
test_points = sample(data_sp$ID[data_sp$model_obs], N, replace=F)
data_sp$partition = ifelse(data_sp$ID %in% test_points, 'test', ifelse(data_sp$model_obs, 'fit', 'excluded'))

# Make partition into correct levels
data_sp$partition = factor(data_sp$partition, levels=c('fit','test','excluded'))
data_ll$partition = data_sp$partition

filename = paste('model_data_',focal_scale,'km_',focal_n,'_records', sep='')

png(paste('./Spatial Model Data/map_',filename,'.png', sep=''), height=1000, width=1000)
spplot(data_sp, 'partition',  main='', panel=function(x,y,subscripts,...){
	panel.pointsplot(x,y,...)
	sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
	SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
		scale = 1000, fill=c("transparent","black"))
	sp.text(c(1500,-2700), "0")
	sp.text(c(2500,-2700), "1000")
	sp.text(c(2000,-2200), "km")
		
	}, cex=0.6, pch=15, col.regions = c('green','blue','grey80'), auto.key=F,
	key=list(x=.1,y=.5, corner=c(0,.5), title='',
		rectangles=list(col=c('green','blue','grey80'), size=1, border='transparent'),
		text=list(paste(levels(data_sp$partition), 'obs')))
)
dev.off()


# Save data frame for analysis
#writeOGR(data_ll, 'Spatial Model Data', filename, 'ESRI Shapefile')
#write.csv(data_ll@data, paste('./Spatial Model Data/', filename, '.csv', sep=''), row.names=F)

### Fit models for Parmeliaceae and Physciaceae separately
library(geoR)
library(spBayes)

## Explore data structure

# Distributions of response and residuals - all appear normal
par(mfrow=c(1,2))
hist(data_ll$S_parm);hist(data_ll$S_phys)

use_vars = names(data_ll)[grep('mean|cv',names(data_ll))]
use_obs = data_ll$partition=='fit'

parm_mod0 = lm(S_parm~., data = data_ll@data[use_obs,c('S_parm',use_vars)])
plot(parm_mod0)
hist(resid(parm_mod0))
phys_mod0 = lm(S_phys~., data = data_ll@data[use_obs,c('S_phys',use_vars)])
plot(phys_mod0)
hist(resid(phys_mod0))

# Convert to geodata object (used by geoR package) - check that column numbers are correct
parm_gd = as.geodata(data_sp[data_sp$partition=='fit',], data.col=4, covar.col=c(6:10, 16:20))
phys_gd = as.geodata(data_sp[data_sp$partition=='fit',], data.col=5, covar.col=c(6:10, 16:20))
plot(parm_gd)
plot(phys_gd)

# Empirical variograms
bin_parm = variog(parm_gd, uvec=seq(0,7200,by=300))
bin_phys = variog(phys_gd, uvec=seq(0,7200,by=300))

par(mfrow = c(1,2))
plot(bin_parm, main = "Parmeliaceae")
plot(bin_phys,  main = "Physciaceae")

# plot variogram models
fit_parm = likfit(parm_gd, ini = c(100,2000), nug = 1)
fit_phys = likfit(phys_gd, ini = c(100,2000), nug = 1)

par(mfrow = c(1,2))
plot(bin_parm, main = "Parmeliaceae")
lines.variomodel(cov.model="exp", cov.pars=c(120,1708), nug=0, max.dist=3000)
lines.variomodel(cov.model="mat", cov.pars=c(120,1708), nug=0, kappa=1, max.dist=3000, lty=2)
lines.variomodel(cov.model="sph", cov.pars=c(120,1708), nug=0, max.dist=3000, lwd=2)
plot(bin_phys,  main = "Physciaceae")
lines.variomodel(cov.model="exp", cov.pars=c(126,2282), nug=0, max.dist=3000)
lines.variomodel(cov.model="mat", cov.pars=c(126,2282), nug=0, kappa=1, max.dist=3000, lty=2)
lines.variomodel(cov.model="sph", cov.pars=c(126,2282), nug=0, max.dist=3000, lwd=2)
lines.variomodel(cov.model="gaus", cov.pars=c(126,2282), nug=0, max.dist=3000, lty=3)
# use spherical model for parm and exponential model for phys


### Bayesian estimates of model using spBayes

# Parmeliaceae model
dat = data.frame(richness=parm_gd$data, parm_gd$covariate)
dat[,2:ncol(dat)] = scale(dat[,2:ncol(dat)], center=T, scale=F)

# Define priors
p = ncol(dat) # num env covariates
phi_range = 3/c(3500, 100) # min and max of phi (range)
sigma_mean = var(dat$richness)
priors =  list('beta.norm'=list(rep(0,p), diag(1000,p)),
	'phi.unif'=phi_range, 'sigma.sq.ig'=c(2, sigma_mean))

# Define starting values and sampler update increments
starting = list('beta'=rep(0,p), 'phi'=3/1000, 'sigma.sq'=sigma_mean)
tuning = list('phi'=.1, 'sigma.sq'=.001)

# Fit model - tau = 0 (no nugget)
parm_mod = spLM(richness~., data = dat, coords = parm_gd$coords,
	starting=starting, priors=priors, tuning=tuning, cov.model='spherical',
	n.samples=5000, n.report=100, verbose=T)

# Check convergence
c1 = spRecover(parm_mod, start=1, thin=10)
c2 = spRecover(parm_mod, start=1, thin=10)
combinedchains = mcmc.list(c1$p.beta.recover.samples, c2$p.beta.recover.samples)
gelman.diag(combinedchains)
gelman.plot(combinedchains) # Use burnin of at least 5000. Convergence after 5000.
combinedchains = mcmc.list(c1$p.theta.recover.samples, c2$p.theta.recover.samples)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

# Physciaceae model
dat = data.frame(richness=phys_gd$data, phys_gd$covariate) 
dat[,2:ncol(dat)] = scale(dat[,2:ncol(dat)], center=T, scale=F)

# Define priors
p = ncol(dat) # num env covariates
sigma_mean = var(dat$richness) # sets mean of inverse gamma prior to empirical variance
phi_range = 3/c(3500, 100) # sets support of uniform prior to effective range: -log(0.05) / min and max of range
priors =  list('beta.norm'=list(rep(0,p), diag(1000,p)),
	'phi.unif'=phi_range, 'sigma.sq.ig'=c(2, sigma_mean))

# Define starting values and sampler update increments
starting = list('beta'=rep(0,p), 'phi'=3/1000, 'sigma.sq'=sigma_mean)
tuning = list('phi'=.1, 'sigma.sq'=.001)

# Fit model - tau = 0 (no nugget)
phys_mod = spLM(richness~., data = dat, coords = phys_gd$coords,
	starting=starting, priors=priors, tuning=tuning, cov.model='exponential',
	n.samples=100000, n.report=10000, verbose=T)

# Recover samples from runs
parm_samps = spRecover(parm_mod, start=1, thin=10)
phys_samps = spRecover(phys_mod, start=1, thin=10)

# Plot traces of samples
plot(parm_samps$p.beta.recover.samples)
plot(parm_samps$p.theta.recover.samples)
plot(phys_samps$p.beta.recover.samples)
plot(phys_samps$p.theta.recover.samples)

# Define burnin samples to discard
burnin = 1000 # this actually starts at 10000 of 100000
end_samp = 100000/10

# Examine autocorrelation to determine whether samples should be thinned
acf(parm_samps$p.theta.recover.samples[burnin:10000,1]) # 18
acf(parm_samps$p.theta.recover.samples[burnin:10000,2]) # 12
acf(parm_samps$p.beta.recover.samples[burnin:10000,1]) # 1
acf(phys_samps$p.theta.recover.samples[burnin:10000,1], lag.max=100) # 42
acf(phys_samps$p.theta.recover.samples[burnin:10000,2], lag.max=100) # 42
acf(phys_samps$p.beta.recover.samples[burnin:10000,1]) # 1
abline(h=0.1)

# Define thinning interval (samples to use)
thin = 20
keep_samps = seq(burnin, end_samp, by=thin)

## Examine data from models run on computing cluster
# Three chains of each model were run on the KillDevil computing cluster using scripts fit_parm_spLM.R and fit_phys_spLM.R

# run this code for each of the 3 chains
# samples have no burnin but have already een thinned by 200
# Note that phys models were incorrectly run with parmeliaceae richness prior to run 7
load('./Spatial Model Data/parm_spLM_run1.RData')
parm_c1 = parm_samps # used when examining convergence from several shorter runs
parm_c2 = parm_samps
parm_c3 = parm_samps

load('./Spatial Model Data/phys_spLM_run1.RData')
#phys_c1 = phys_samps
#phys_c2 = phys_samps
#phys_c3 = phys_samps

# Check convergence among chains 
c1 = phys_c1 #parm_c1
c2 = phys_c2 #parm_c2
c3 = phys_c3 #parm_c3

beta_combo = mcmc.list(c1$p.beta.recover.samples, c2$p.beta.recover.samples, c3$p.beta.recover.samples)
theta_combo = mcmc.list(c1$p.theta.recover.samples, c2$p.theta.recover.samples, c3$p.theta.recover.samples)

gelman.plot(beta_combo, mfrow=c(2,5)) # all less than 1.1 after i=200
gelman.plot(theta_combo, mfrow=c(2,5))  # Use burnin of at least 5000. Convergence after 5000.
# all parameters converge with shrink factor < 1.1 after i=200 for Parmeliaceae and i=300 for Physciaceae
# which corresponds to a burnin of 40000 or 60000.

# Plot traces and check autocorrelation
par(mfrow=c(1,2))
plot(beta_combo, auto.layout=F, ask=T)
par(mfrow=c(1,2))
plot(theta_combo, auto.layout=F, ask=T)

burnin=301
last_samp=500

for(i in 1:11){
	acf(c1$p.beta.recover.samples[burnin:last_samp,i], main=colnames(c1$p.beta.recover.samples)[i])
	invisible(readline(prompt="Press [enter] to continue"))
}
acf(c1$p.theta.recover.samples[burnin:last_samp,1])
acf(c1$p.theta.recover.samples[burnin:last_samp,2])
# No autocorrelation between samples. Ok to use every sample.

# Extract samples to be used for parameter estimation into single mcmc object
parm_samps3 = cbind(beta_combo[[1]][burnin:last_samp,], theta_combo[[1]][burnin:last_samp,])
for(i in 2:3){
	new_dat = cbind(beta_combo[[i]][burnin:last_samp,], theta_combo[[i]][burnin:last_samp,])
	parm_samps3 = rbind(parm_samps3, new_dat)
}

phys_samps3 = cbind(beta_combo[[1]][burnin:last_samp,], theta_combo[[1]][burnin:last_samp,])
for(i in 2:3){
	new_dat = cbind(beta_combo[[i]][burnin:last_samp,], theta_combo[[i]][burnin:last_samp,])
	phys_samps3 = rbind(phys_samps3, new_dat)
}

# Examine correlations between parameters across samples
pairs(parm_samps3)
pairs(phys_samps3)

# Focus on pair of parameters that appear to be correlated
use_samps = phys_samps3
xnum = 12
ynum = 13

plot(use_samps[,xnum], use_samps[,ynum], xlab=colnames(use_samps)[xnum], ylab=colnames(use_samps)[ynum])
abline(h=0,v=0)
xest = quantile(use_samps[,xnum], probs=c(0.025, .5, 0.975))
yest = quantile(use_samps[,ynum], probs=c(0.025, .5, 0.975))
usr = par('usr')
lines(xest, rep(usr[3],3), col=2, lwd=4)
lines(rep(usr[1],3), yest , col=2, lwd=4)

## Re-ran models with 10^6 iterations to be able to use spPredict

# Load data from long runs
load('./Spatial Model Data/parm_spLM_run7.RData')
load('./Spatial Model Data/phys_spLM_run7.RData')

# Combine beta and theta samples so that correlations can be examine
parm_allsamps = cbind(parm_samps$p.beta.recover.samples, parm_samps$p.theta.recover.samples)
phys_allsamps = cbind(phys_samps$p.beta.recover.samples, phys_samps$p.theta.recover.samples)

pairs(parm_allsamps[300:5000,])
pairs(phys_allsamps[300:5000,])

# Define samples to use
thin=200
burnin=300*thin+1

# Create prediction data
parm_gd_test = as.geodata(data_sp[data_sp$partition=='test',], data.col=4, covar.col=c(6:10, 16:20))
phys_gd_test = as.geodata(data_sp[data_sp$partition=='test',], data.col=5, covar.col=c(6:10, 16:20))

# Make prediction
parm_dat = parm_gd_test$covariate
parm_dat = scale(parm_dat, center=T, scale=F)
parm_dat = cbind(rep(1,nrow(parm_dat)), parm_dat) 

phys_dat = phys_gd_test$covariate
phys_dat = scale(phys_dat, center=T, scale=F)
phys_dat = cbind(rep(1,nrow(phys_dat)), phys_dat) 

parm2parm = spPredict(parm_mod, parm_gd_test$coords, parm_dat, start=burnin, thin=thin)
phys2phys = spPredict(phys_mod, phys_gd_test$coords, phys_dat, start=burnin, thin=thin)
parm2phys = spPredict(parm_mod, phys_gd_test$coords, phys_dat, start=burnin, thin=thin)
phys2parm = spPredict(phys_mod, parm_gd_test$coords, parm_dat, start=burnin, thin=thin)

## Calculate probability of observed richness under posterior distributions
calc_p = function(x, dist){
	sapply(1:length(x), function(i){
		greater = sum(dist[i,] >= x[i])
		lesser = sum(dist[i,] <= x[i])
		2*min(greater, lesser)/ncol(dist)	
	})
}

P_parm2parm = calc_p(parm_gd_test$data, parm2parm$p.y.predictive.samples)
P_phys2phys = calc_p(phys_gd_test$data, phys2phys$p.y.predictive.samples)
P_parm2phys = calc_p(phys_gd_test$data, parm2phys$p.y.predictive.samples)
P_phys2parm = calc_p(parm_gd_test$data, phys2parm$p.y.predictive.samples)

# Plot histograms of P-values for each model-data combination
par(mfrow=c(2,2))
hist(P_parm2parm,xlim=c(0,1))
hist(P_parm2phys,xlim=c(0,1))
hist(P_phys2parm,xlim=c(0,1))
hist(P_phys2phys,xlim=c(0,1))


## Calculate residual deviation of observed richness from median of predicted posterior
calc_res = function(x, dist){
	meds = apply(dist, 1, median)
	meds - x	
}


res_parm2parm = calc_res(parm_gd_test$data, parm2parm$p.y.predictive.samples)
res_phys2phys = calc_res(phys_gd_test$data, phys2phys$p.y.predictive.samples)
res_parm2phys = calc_res(phys_gd_test$data, parm2phys$p.y.predictive.samples)
res_phys2parm = calc_res(parm_gd_test$data, phys2parm$p.y.predictive.samples)


## Map residual deviation and probabilities

# Make into spatial data
modcomp_df = data.frame(P_parm2parm, res_parm2parm, P_phys2phys, res_phys2phys, P_parm2phys, res_parm2phys, P_phys2parm, res_phys2parm)
rownames(modcomp_df) = rownames(parm_gd_test$coords)
modcomp_sp = SpatialPointsDataFrame(parm_gd_test$coords, modcomp_df, proj4string=CRS( proj4string(data_sp)))


colcuts = c(0,0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1)
ncuts=length(colcuts) - 1
mycolramp = colorRampPalette(mycol)(ncuts)

pdf('spatial model comparison P.pdf', height=8, width=8)
for(var in c('P_parm2parm', 'P_phys2phys', 'P_parm2phys','P_phys2parm')){
print(
spplot(modcomp_sp, var,  main=var, panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(space='right', corner=c(0,.5), title='P',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(as.character(colcuts[2:length(colcuts)]))
	)
)
)
}
dev.off()

range(c(res_parm2parm, res_phys2phys, res_parm2phys,res_phys2parm))
colcuts = seq(-50, 50, 10)
ncuts=length(colcuts) - 1
mycolramp = colorRampPalette(mycol)(ncuts)

pdf('spatial model comparison residuals.pdf', height=8, width=8)
for(var in c('res_parm2parm', 'res_phys2phys', 'res_parm2phys','res_phys2parm')){
print(
spplot(modcomp_sp, var,  main=var, panel=function(x,y,subscripts,...){
		panel.pointsplot(x,y,...)
		sp.polygons(nam_laea, fill='transparent', col="#4c4c4ccc")
		
		SpatialPolygonsRescale(layout.scale.bar(), offset = c(1500,-2500), 
			scale = 1000, fill=c("transparent","black"))
		sp.text(c(1500,-2700), "0")
		sp.text(c(2500,-2700), "1000")
		sp.text(c(2000,-2200), "km")
	}, cuts=colcuts, cex=0.6, pch=15, col.regions = mycolramp, auto.key=F,
	key=list(space='right', corner=c(0,.5), title='# Species',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(as.character(colcuts[2:length(colcuts)]))
	)
)
)
}
dev.off()

# Plot observed richness versus predicted distributions

plot_obs_pred = function(x, dist){
	confints = apply(dist, 1, function(x) quantile(x, probs=c(0.025, 0.975)))
	meds = apply(dist, 1, median)
	myrange = range(confints)
	plot(x, meds, pch=16, xlim=myrange, ylim=myrange, xlab='Observed richness', ylab='Predicted richness')
	for(i in 1:length(x)){
		lines(rep(x[i],2), confints[,i])
	}
	abline(0,1, col=2)
}

pdf('observed and predicted richness in parm and phys models.pdf', height=8, width=8)
par(mfrow=c(2,2))
par(mar=c(4,6,6,1))
plot_obs_pred(parm_gd_test$data, parm2parm$p.y.predictive.samples)
mtext('Parmeliaceae Model', 3, 1)
mtext('Parmeliaceae Richness', 2, 5)
plot_obs_pred(parm_gd_test$data, phys2parm$p.y.predictive.samples)
mtext('Physciaceae Model', 3, 1)
par(mar=c(6,6,4,1))
plot_obs_pred(phys_gd_test$data, parm2phys$p.y.predictive.samples)
mtext('Physciaceae Richness', 2, 5)
plot_obs_pred(phys_gd_test$data, phys2phys$p.y.predictive.samples)
dev.off()

## Calculate likelihood of observed data under both models using predicted distributions
omit0 = which(P_parm2parm==0)
prod(P_parm2parm[-omit0])

## Try calculating R2?
SS = sum(res_parm2parm^2)


## Compare parameter estimates
burnin = 301
lastsamp = 5000

# Naive parameter estimates and credibility intervals marginalizing across all other parameters (assumes independences of MCMC samples across parameters)

parm_ests = apply(parm_allsamps[burnin:lastsamp,], 2, median)
phys_ests = apply(phys_allsamps[burnin:lastsamp,], 2, median)

parm_ints = apply(parm_allsamps[burnin:lastsamp,], 2, function(x) quantile(x, probs=c(0.025, 0.975)))
phys_ints = apply(phys_allsamps[burnin:lastsamp,], 2, function(x) quantile(x, probs=c(0.025, 0.975)))

parm_sig = apply(parm_ints, 2, function(x) prod(x)>0)
phys_sig = apply(phys_ints, 2, function(x) prod(x)>0)

data.frame(parm_ests, parm_sig, t(parm_ints), phys_ests, phys_sig, t(phys_ints))






