library(TESS)

# Change this to the directory containing the dryad package
dir <- '~/repos/mass-extinction-simulation/dryad_package/simulation_5'

# Read in the conifer tree
conifer_tree <- read.nexus(paste(dir,'/data/conifer_40my.tre',sep=''))

# Read in output from one of the empirical analyses
output <- tess.process.output(dir = paste(dir,'/empirical_analysis_output',sep=''), tree = conifer_tree)

# Simulation settings
N <- length(conifer_tree$tip.label)
time <- max(branching.times(conifer_tree))
t_1 <- 250
t_2 <- c(200,150,100,50)
survival_probability <- c(0.1,0.1)
trees <- 1:100

grid <- expand.grid(t_2=t_2,tree=trees)

invisible(apply(grid,1,function(row){
  
  second_event_time <- row[[1]]
  T <- row[[2]]
  
  tree_dir <- paste(dir,'/data/n_',second_event_time,'/tree_',T,sep='')
  dir.create(tree_dir,showWarnings=FALSE,recursive=TRUE)
  
  lambda <- sample(output$spec,size=1)
  mu <- sample(output$ext,size=1)
  
  while(mu > lambda){
    lambda <- sample(output$spec,size=1)
    mu <- sample(output$ext,size=1)
  }
  
  tree <- tess.sim.taxa.age(n=1,nTaxa=n,age=t,lambda=lambda,mu=mu,massExtinctionTimes=t-c(t_1,time),massExtinctionSurvivalProbabilities=survival_probability)[[1]]
  save(tree,file=paste(tree_dir,'/tree.Robj',sep=''))
  
  pdf(file=paste(tree_dir,"/tree_figure.pdf",sep=""),height=6,width=12)
  par(mfrow=c(1,2))
  plot(tree,show.tip.label=FALSE,edge.width=1,x.lim=c(0,t))
  abline(v=t-c(t_1,second_event_time),col="purple",lwd=1)
  matplot(x=c(-t,0),y=cbind(rep(lambda,2),rep(mu,2)),type="l",lty=1,col=c("blue","red"),xlab="time",ylab="rate",lwd=1,lend=2,ylim=c(0,max(cbind(lambda,mu))))
  abline(v=-c(t_1,second_event_time),col="purple",lwd=1,lend=2)
  axis(3,at=-c(t_1,second_event_time),labels=round(survival_probability,5),las=2,tick=FALSE,line=-0.5)
  dev.off()
  
}))

