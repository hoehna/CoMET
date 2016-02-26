#########################
## Run a tess analysis ##
#########################

## Argument 1 is the working directory.
## Argument 2 is the number of cores for mclapply to use.
## Argument 3 is the job ID.
## Argument 4 is the tree size.
## Argument 5 is whether to use empirical hyperprior.
## Argument 6 is the speciation rate prior mean (in log space).
## Argument 7 is the speciation rate prior standard deviation (in log space).
## Argument 8 is the extinction rate prior mean (in log space).
## Argument 9 is the extinction rate prior standard deviation (in log space).
## Argument 10 is the diversification-rate-shift rate.
## Argument 11 is the mass-extinction rate.

## Load packages.
library(TESS)
library(parallel)

## Get the arguments
args <- commandArgs(TRUE)

# cat(args)

dir <- args[1]
MC <- as.numeric(args[2])
ID <- as.numeric(args[3])
N <- as.numeric(args[4])
EMP <- as.logical(args[5])
LM <- as.numeric(args[6])
LS <- as.numeric(args[7])
MM <- as.numeric(args[8])
MS <- as.numeric(args[9])
S <- eval(parse(text=args[10]))
M <- eval(parse(text=args[11]))

invisible(mclapply(1:100,function(x){

  ## Get the tree directory and load the tree
  tree_dir <- paste(dir,'/data/n_',N,'/tree_',x,'/tree.Robj',sep='')
  load(tree_dir)
  
  ## Get the output directory
  out_dir <- paste(dir,'/output/job_',ID,'/tree_',x,sep='')
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  
  ## Perform the analysis
  smallest_ESS <- 0
  number_of_tries <- 1
  
  while ( smallest_ESS < 250 ) { # Redo any analyses with ESS < 250
    
    if ( number_of_tries > 5 ) { # If more than 5 retries, something is seriously problematic. Check the tree and MCMC output!
      
      cat('Too many retries. Check tree and MCMC output manually!',file=paste(out_dir,'/warning.txt',sep=''))
      
    } else {
      
      tess.analysis(tree = tree,
                    empiricalHyperPriors = EMP,
                    initialSpeciationRate = 1,
                    speciationRatePriorMean = LM,
                    speciationRatePriorStDev = LS,
                    initialExtinctionRate = 0.5,
                    extinctionRatePriorMean = MM,
                    extinctionRatePriorStDev = MS,
                    estimateNumberRateChanges = TRUE,
                    numExpectedRateChanges = S,
                    pMassExtinctionPriorShape1 = 2,
                    pMassExtinctionPriorShape2 = 18,
                    estimateMassExtinctionTimes = TRUE,
                    numExpectedMassExtinctions = M,
                    MAX_ITERATIONS = 1e6,
                    THINNING = 100,
                    MIN_ESS = 500,
                    CONVERGENCE_FREQUENCY = 1000,
                    dir = out_dir)
      
      output <- read.table(paste(out_dir,'/samples_numCategories.txt',sep=''),header=TRUE)
      smallest_ESS <- min(effectiveSize(output[,3:5]))
      number_of_tries <- number_of_tries + 1
      
    }
    
  }
  
},mc.cores=MC,mc.preschedule=FALSE))

# dir <- '~/repos/mass-extinction-simulation/dryad_package/simulation_1'
# MC <- 16
# ID <- 1
# N <- 100
# EMP <- FALSE
# LM <- -0.09116078
# LS <- 0.4269913
# MM <- -3.824846
# MS <- 1.744856
# S <- 0.1
# M <- 0.1
# 
# invisible(lapply(1:100,function(x){
# 
#   ## Get the tree directory and load the tree
#   tree_dir <- paste(dir,'/data/n_',N,'/tree_',x,'/tree.Robj',sep='')
#   load(tree_dir)
# 
#   ## Get the output directory
#   out_dir <- paste(dir,'/output/job_',ID,'/tree_',x,sep='')
#   dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
#   
#   ## Perform the analysis
#   smallest_ESS <- 0
#   number_of_tries <- 1
#   
#   while ( smallest_ESS < 250 ) { # Redo any analyses with ESS < 250
#     
#     if ( number_of_tries > 5 ) { # If more than 5 retries, something is seriously problematic. Check the tree and MCMC output!
#       
#       cat('Too many retries. Check tree and MCMC output manually!',file=paste(out_dir,'/warning.txt',sep=''))
#       
#     } else {
#       
#       tess.analysis(tree = tree,
#                     empiricalHyperPriors = EMP,
#                     initialSpeciationRate = 1,
#                     speciationRatePriorMean = LM,
#                     speciationRatePriorStDev = LS,
#                     initialExtinctionRate = 0.5,
#                     extinctionRatePriorMean = MM,
#                     extinctionRatePriorStDev = MS,
#                     estimateNumberRateChanges = TRUE,
#                     numExpectedRateChanges = S,
#                     pMassExtinctionPriorShape1 = 2,
#                     pMassExtinctionPriorShape2 = 18,
#                     estimateMassExtinctionTimes = TRUE,
#                     numExpectedMassExtinctions = M,
#                     MAX_ITERATIONS = 1e6,
#                     THINNING = 100,
#                     MIN_ESS = 500,
#                     CONVERGENCE_FREQUENCY = 1000,
#                     dir = out_dir)
# 
#       output <- read.table(paste(out_dir,'/samples_numCategories.txt',sep=''),header=TRUE)
#       smallest_ESS <- min(effectiveSize(output[,3:5]))
#       number_of_tries <- number_of_tries + 1
#       
#     }
#     
#   }
# 
# }))



















