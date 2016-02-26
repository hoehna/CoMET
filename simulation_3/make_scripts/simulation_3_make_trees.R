library(TESS)

# Change this to the directory containing the dryad package
dir <- '~/repos/mass-extinction-simulation/dryad_package/simulation_3'

# Simulation settings
time <- 10
lambda_mean <- 1
lambda_sd <- exp(0.2)
mu_mean <- 0.5
mu_sd <- exp(0.2)
diversification_shift_rate <- 0
mass_extinction_rate <- 1
mass_extinction_alpha <- 2
mass_extinction_beta <- 18

# Parameter combinations
N <- c(100,200,400,800)
trees <- 1:100

grid <- expand.grid(N=N,tree=trees)

invisible(apply(grid,1,function(row){
  
  N <- row[[1]]
  T <- row[[2]]
  
  tree_dir <- paste(dir,'/data/n_',N,'/tree_',T,sep='')
  dir.create(tree_dir,showWarnings=FALSE,recursive=TRUE)
  
  # Make speciation rate shifts.
  num_lambda_shifts = rpois(n=1,lambda=diversification_shift_rate)
  times_lambda_shifts = sort(runif(n=num_lambda_shifts,min=0,max=time))
  rates_lambda_shifts = rlnorm(n=num_lambda_shifts+1,meanlog=log(lambda_mean),sdlog=log(lambda_sd))
  lambda_function = function(x) rates_lambda_shifts[1+findInterval(x,times_lambda_shifts)]
  
  # Make extinction rate shifts.
  num_mu_shifts = rpois(n=1,lambda=diversification_shift_rate)
  times_mu_shifts = sort(runif(n=num_mu_shifts,min=0,max=time))
  rates_mu_shifts = rlnorm(n=num_mu_shifts+1,meanlog=log(mu_mean),sdlog=log(mu_sd))
  mu_function = function(x) rates_mu_shifts[1+findInterval(x,times_mu_shifts)]
  
  # Make mass extinction events.
  num_mass = rpois(n=1,lambda=mass_extinction_rate)
  times_mass = sort(runif(n=num_mass,min=0,max=time))
  survival_probs = rbeta(n=num_mass,shape1=mass_extinction_alpha,shape2=mass_extinction_beta)
  
  # Simulate the tree, store diversification process parameters in the tree
  tree <- tess.sim.taxa.age(n=1,nTaxa=N,age=time,lambda=lambda_function,mu=mu_function,
                            massExtinctionTimes=times_mass,massExtinctionSurvivalProbabilities=survival_probs)[[1]]
  
  event_list <- list(times_lambda_shifts = times_lambda_shifts,
                     rates_lambda_shifts = rates_lambda_shifts,
                     times_mu_shifts = times_mu_shifts,
                     rates_mi_shifts = rates_mu_shifts,
                     mass_extinction_times = times_mass,
                     survival_probs = survival_probs)

  tree$event_list <- event_list

  # Save the tree
  save(tree,file=paste(tree_dir,'/tree.Robj',sep=''))
  
  # Make a figure of the tree and the process
  times <- seq(0,time,length.out=10000)
  lambdas <- lambda_function(times)
  mus <- mu_function(times)
  
  pdf(file=paste(path,"/tree_figure.pdf",sep=""),height=6,width=12)
  par(mfrow=c(1,2))
  plot(tree,show.tip.label=FALSE,edge.width=1,x.lim=c(0,time))
  abline(v=times_lambda_shifts,col="purple",lwd=1)
  abline(v=times_mu_shifts,col="red",lwd=1)
  abline(v=times_mass,col="green",lwd=1)
  matplot(x=times,y=cbind(lambdas,mus),type="l",lty=1,col=c("purple","red"),xlab="time",ylab="rate",lwd=1,lend=2,ylim=c(0,max(cbind(lambdas,mus))))
  abline(v=times_mass,col="green",lwd=1,lend=2)
  axis(3,at=times_mass,labels=round(survival_probs,5),las=2,tick=FALSE,line=-0.5)
  dev.off()
  
}))

