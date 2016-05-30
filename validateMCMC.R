# check MCMC run for convergence and create report of plots
library(coda)
runsPath <- "~/data/IES/analysis/asr/"
runDirs <- c("run1/", "run2/")
runs <- paste0(runsPath, runDirs)
mcmcL <- list()
counter <- 1
for(i in runs){
  #  i <- runs[1]
  traceI <- read.table(paste0(i,"parameters1.log"), header = TRUE)
  thin <- traceI$Iteration[2] - traceI$Iteration[1]
  mcmcL[[counter]] <- mcmc(data = subset(traceI, select = -Iteration), thin = thin)
  counter <- counter + 1
}


## The parameter space was broadly sampled

plot(mcmcL[[1]])
plot(mcmcL[[2]])

## The time series had minimal autocorrelation

autocorr.plot(mcmcL[[1]])
autocorr.plot(mcmcL[[2]])

## The runs converged

gelman.plot(mcmcL)

