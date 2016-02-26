# Data and priors used for the analysis
N <- c(100,200,400,800)
extinction_rate <- c(0.1,0.5,0.9,'empirical')
diversification_shift_rate <- c(0.1,'log(2)',2,5)
mass_extinction_rate <- c(0.1,'log(2)',2,5)

# Change this to the number of processors mclapply will use
num_cores <- 16 

# Change this to the directory containing the dryad package
dir <- '~/repos/mass-extinction-simulation/dryad_package/simulation_2' 

# This grid contrains all combinations of data and prior settings
grid <- expand.grid(N=N,
                    extinction_rate=extinction_rate,
                    diversification_shift_rate=diversification_shift_rate,
                    mass_extinction_rate=mass_extinction_rate)
grid$ID <- 1:nrow(grid)

# For each combination of data set and prior settings, write a script to analyze all 100 trees.
invisible(apply(grid,1,function(row){

  N <- row[[1]]
  E <- row[[2]]
  S <- row[[3]]
  M <- row[[4]]

  ## Get the empirical hyperprior bool
  if ( E == 'empirical' ) {
    empirical_hyperprior <- TRUE
    E <- 0.5
  } else {
    empirical_hyperprior <- FALSE
    E <- as.numeric(E)
  }

  ## Compute the parameters of the lognormal distribution on speciation and extinction rates
  m_lambda = 1
  v_lambda = 0.2
  mu_lambda = log((m_lambda^2)/sqrt(v_lambda+m_lambda^2))
  var_lambda = log(1+v_lambda/(m_lambda^2))

  m_mu = m_lambda*E
  v_mu = 0.2
  mu_mu = log((m_mu^2)/sqrt(v_mu+m_mu^2))
  var_mu = log(1+v_mu/(m_mu^2))

  ## Get the path to the script file, and make the script file.
  script_file <- paste(dir,'/run_scripts/job_',as.numeric(row[[5]]),'.sh',sep='')
  file.create(script_file,showWarnings=FALSE,recursive=TRUE)

  ## Now, write the script file.
  cat('Rscript ',dir,'/run_scripts/runTessAnalysis.R ',
      dir,' ',                      # Argument 1
      num_cores,' ',                # Argument 2
      as.numeric(row[[5]]),' ',     # Argument 3
      N,' ',                        # Argument 4
      empirical_hyperprior,' ',     # Argument 5
      mu_lambda,' ',                # Argument 6
      sqrt(var_lambda),' ',         # Argument 7
      mu_mu,' ',                    # Argument 8
      sqrt(var_mu),' ',             # Argument 9
      S,' ',                        # Argument 10
      M,                            # Argument 11
      sep='',append=TRUE,file=script_file)

}))






