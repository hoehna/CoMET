# Data and priors used for the analysis
extinction_rate <- c('empirical')
diversification_shift_rate <- c(0.1,log(2),2,5)
mass_extinction_rate <- c(0.1,log(2),2,5)
runs <- 1:4

# Change this to the number of processors mclapply will use
num_cores <- 4

# Change this to the directory containing the dryad package
dir <- '~/repos/mass-extinction-simulation/dryad_package/conifers_empirical_analysis'

# Load the data
tree <- read.nexus(paste(dir,'/data/conifer_40my.tre',sep=''))
rho <- (tree$Nnode+1)/630

# This grid contrains all combinations of data and prior settings
grid <- expand.grid(runs=runs,
                    extinction_rate=extinction_rate,
                    diversification_shift_rate=diversification_shift_rate,
                    mass_extinction_rate=mass_extinction_rate)
grid$ID <- 1:nrow(grid)

# For each combination of prior settings, analyze the confier data
invisible(mclapply(1:nrow(grid),function(x){
  
  row <- grid[x,]
  
  S <- row[[3]]
  M <- row[[4]]
  R <- row[[1]]
  
  this_dir <- paste(dir,'/output/shiftrate_',S,'/merate_',M,'/run_',R,sep='')
  dir.create(this_dir,showWarnings=FALSE,recursive=TRUE)
  
  tess.analysis(tree,initialSpeciationRate=1,initialExtinctionRate=0.5,
                empiricalHyperPriors=TRUE,samplingProbability=rho,
                numExpectedRateChanges = S,pMassExtinctionPriorShape1 = 2.5,
                pMassExtinctionPriorShape2 = 7.5,numExpectedMassExtinctions = M,
                MAX_ITERATIONS = 2e6,MIN_ESS = 500,dir = this_dir)
    
},mc.cores=num_cores))