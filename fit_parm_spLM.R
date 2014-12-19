## This script estimates parameters for a spatial model of regional parmeliaceae richness

### Read in Data and Packages

library(spBayes, lib.loc='~/Rlibs')
library(geoR)
library(sp)
library(rgdal)

plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

data_ll = readOGR('Spatial Model Data','model_data_250km_100_records')
data_sp = spTransform(data_ll, CRS(plot_prj))


# Get arguments passed from R CMD BATCH
args = commandArgs(trailingOnly = F)
myarg = args[length(args)] # my argument is at the end
run = as.numeric(sub("-","",myarg))
print(run)


### Code

# Define variables and observations to use
use_vars = names(data_ll)[grep('mean|cv',names(data_ll))]
use_obs = data_ll$partition=='fit'

# Convert to geodata object (used by geoR package) - check that column numbers are correct
parm_gd = as.geodata(data_sp[data_sp$partition=='fit',], data.col=4, covar.col=c(6:10, 16:20))

# Make data frame
dat = data.frame(richness=parm_gd$data, parm_gd$covariate)
dat[,2:ncol(dat)] = scale(dat[,2:ncol(dat)], center=T, scale=F)

## Define priors
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
	n.samples=1000000, n.report=10000, verbose=T)

# Recover samples from runs
parm_samps = spRecover(parm_mod, start=1, thin=200)

# Save model
save(parm_mod, parm_samps, file=paste('parm_spLM_run', run, '.RData', sep=''))







