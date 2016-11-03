library(coda)
library(mcmcplots)
library(rjags)
library(R6)


JAGS <- R6Class("JAGS",
  public = list(
    modelString = NA,
    params = NA,
    config = list(n.chains=3,n.adapt=2400,n.iter.update=2400,n.iter.samples=2400, thin=3,seed=1111),
    mcmc = NA,
    l = list(),
    
    initialize = function(modelString=NA,params=NA,n.chains=3,config=NA,l=list()) {
      self$modelString <- modelString
      self$params <- params
      if (any(!is.na(config))) {
        for (item in config) {self$config[[item]] <- config[[item]]  }
      }
      self$add(l)
    },
    
    add = function(l) {
      self$l <- c(self$l, l)
    },
    
    run = function(data, period=NA, getNIter=TRUE, runDiag=TRUE, showMcmcPlot=TRUE, showGelmanPlot=FALSE) {
      model <- jags.model(textConnection(self$modelString), data=data, n.adapt=self$config$n.adapt, n.chains=self$config$n.chains )
      set.seed(self$seed)
      update(model, n.iter=self$config$n.iter.update)
      self$mcmc <- coda.samples(model=model, variable.names=self$params, n.iter=self$config$n.iter.samples, thin=self$config$thin)
      print(summary(self$mcmc))
      if (getNIter)
        self$config$n.iter <- self$getNIterRafDiag(self$mcmc) 
      if (runDiag) 
        self$getHeidConvDiag(self$mcmc)
      if (showMcmcPlot) 
        for (param in self$params)  {
          cat(sprintf('plotting param %s at period %s ',param, period))
          self$mcmcPlot(self$mcmc, param, period)
        }
      if (showGelmanPlot)
        gelman.plot(self$mcmc)
      return(self$mcmc)
    },
    
    getNIterRafDiag = function(mcmc.output, q=0.025, r=0.005, s=0.95) {
      raf.test <- try(raftery.diag(c(mcmc.output), q=q, r=r, s=s), silent=T)
      if (inherits(raf.test, 'try-error'))
        stop(sprintf('the error level r=%s is too high. Try again with smaller r.',r))
      nmin <- as.numeric(raf.test$resmatrix[2])
      cat(sprintf('Raftery & Lewis N_min interations = %s',nmin))
      return(nmin)
    },
    
    getHeidConvDiag = function(mcmc.output, eps=0.1, pvalue=0.05) {
      cat('\nHeidelberg-Welch diagnostic:\n')
      h <- heidel.diag(mcmc.output, eps=eps, pvalue=pvalue)
      print(h)
      return(all(sapply(h, function(x)x[1,c('stest','htest')]==1)))
    },
    
    mcmcPlot = function(mcmc.output, param=NA, period=NA) {
      n.chains <- length(mcmc.output)
      par(mfrow=c(2,2),mar=c(4,3,2,2))
      mcmcplots::denplot(mcmc.output,style = 'plain',auto.layout = F) #,main=sprintf("Density of %s (t=%s)",param,period))
      mcmcplots::traplot(mcmc.output,style = 'plain',auto.layout = F) #,main=sprintf("Trace of %s (t=%s)",param,period))
      mcmcplots::rmeanplot(mcmc.output,style = 'plain',auto.layout = F) #,main=sprintf("Thinned Running Mean of %s (t=%s)",param,period))
      mcmcplots::autplot1(mcmc.output, chain=n.chains,style = 'plain') #, main=sprintf("Autocorrelation of %s (t=%s)",param,period))
    }
  )
)


###
## EXAMPLE
###
# getModelStr <- function() {
#   "model{
#     for (i in 1:n) {
#       z[i] ~ dbern(q)
#       th1[i] <- p2*(v1 + omega*sig1*z[i])
#       th2[i] <- p1*(v2 + omega*sig2*z[i])
#       s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
#       G[i] ~ dbinom( s[i], L )
#     }
#     epsilon ~ dgamma(shapet,ratet)
#     q ~ dbeta(h1t,h2t)
#   }"
# }
# 
# 
# jg <- JAGS$new(getModelStr(), c('q','epsilon'), list(name='Test1',pds=2))
# . <- jg$run(data.list, period=1)







