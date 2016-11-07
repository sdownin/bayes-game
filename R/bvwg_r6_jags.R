library(coda)
library(mcmcplots)
library(RColorBrewer)
library(MASS)
library(rjags)
library(R2jags)
library(R6)
library(psych)

#######################
## EXAMPLE
#######################
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
# jg <- JAGS$new(getModelStr(), c('q','epsilon'), list(name='Test1',pds=2))
# jg$run(data.list, period=1 )
# jg$bivarPlot('q','epsilon')
######################

JAGS <- R6Class(
  "JAGS",
  public = list(
    probs = c(.005,.025,.5,.975,.995),
    modelstring = NA,
    params = NA,
    config = list(n.chains=3,n.adapt=2400,n.iter.update=2400,n.iter.samples=2400, thin=3,seed=1111),
    mcmc = list(),
    est = list(),
    
    initialize = function(modelstring=NA,params=NA,n.chains=3,...) {
      self$modelstring <- modelstring
      self$params <- params
      self$config$n.chains <- n.chains
      items <- list(...)
      for (item in names(items)){
        if (item %in% names(self$config)) 
          self$config[[item]] <- items[[item]]
        else
          self$addConfig(items[item])
      }
    },
    
    addConfig= function(l) {
      self$config <- c(self$config, l)
    },
    
    run.old = function(data, parallel=T, period=NA, returnOutput=F, getNIter=T, runDiag=T, showMcmcPlot=T, showGelmanPlot=F) {
      model <- jags.model(textConnection(self$modelstring), data=data, n.adapt=self$config$n.adapt, n.chains=self$config$n.chains )
      set.seed(self$seed)
      update(model, n.iter=self$config$n.iter.update)
      self$mcmc <- coda.samples(model=model, variable.names=self$params, n.iter=self$config$n.iter.samples, thin=self$config$thin)
      print(summary(self$mcmc))
      self$setEst()
      
      if (getNIter)
        self$config$n.iter <- self$getNIterRafDiag(self$mcmc) 
      if (runDiag) 
        self$getHeidConvDiag(self$mcmc)
      if (showMcmcPlot) {
        for (param in self$params)  {
          cat(sprintf('plotting param %s at period %s ',param, period))
          self$mcmcPlot(self$mcmc, param, period)
        }
      }
      if (showGelmanPlot)
        gelman.plot(self$mcmc)
      if (returnOutput)
        return(self$mcmc)
    },
    
    ##------------
    run =  function(data, parallel=T, envir=NA, period=NA, returnOutput=F, getNIter=T, runDiag=T, showMcmcPlot=T, showGelmanPlot=F) {
      if (parallel)
        self$mcmc <- self$.runMcmcParallel(data, envir)
      else
        self$mcmc <- self$.runMcmc(data)
      print(summary(self$mcmc))
      self$setEst()
      if (getNIter)
        self$config$n.iter <- self$getNIterRafDiag(self$mcmc) 
      if (runDiag) 
        self$getHeidConvDiag(self$mcmc)
      if (showMcmcPlot) {
        for (param in self$params)  {
          cat(sprintf('plotting param %s at period %s ',param, period))
          self$mcmcPlot(self$mcmc, param, period)
        }
      }
      if (showGelmanPlot)
        gelman.plot(self$mcmc)
      if (returnOutput)
        return(self$mcmc)
    }, 
    
    .runMcmc = function(data)  {
      model <- jags.model(textConnection(self$modelstring), data=data, n.adapt=self$config$n.adapt, n.chains=self$config$n.chains )
      set.seed(self$seed)
      update(model, n.iter=self$config$n.iter.update)
      mcmc.output <- coda.samples(model=model, variable.names=self$params, n.iter=self$config$n.iter.samples, thin=self$config$thin)
      return(mcmc.output)
    }, 
    
    .runMcmcParallel = function(data, envir)   {
      model.file <- "bvwg_r6_jags.bug"
      capture.output(cat(stringr::str_replace_all(self$modelString,"\n","")), file=model.file)
      mcmc.output <- as.mcmc(do.call(jags.parallel,
                                     list(data = data,  parameters.to.save = self$params,
                                          n.chains=4, n.iter=self$config$n.iter.samples, 
                                          n.burnin=ceiling(self$config$n.iter.samples*.3),
                                          model.file=model.file, envir=envir
                                     ), envir = envir
      ))
      return(mcmc.output)
    },
    ##-------------
    
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
    },
    
    bivarPlot = function(varx, vary, chain=NA, cor.title=TRUE, padx=0, pady=.1) {
      if (length(self$mcmc) > 0) {
        if(is.na(chain) | chain=='all') {
          x <- as.vector(unlist(sapply(self$mcmc,function(x)x[, varx])))
          y <- as.vector(unlist(sapply(self$mcmc,function(x)x[, vary])))
        } else {
          x <- as.vector(self$mcmc[[chain]][, varx])
          y <- as.vector(self$mcmc[[chain]][, vary])          
        }
        # density config
        kde <- kde2d(x, y, n=50)
        k <- 10
        my.cols <- rev(brewer.pal(k, "RdYlBu"))
        #
        main <- ''
        if(cor.title) {
          tmp <- cor.test(x,y)[c('estimate','p.value')]
          if(tmp$p.value < 0.001) 
            main <- sprintf('Corr: %.3f (p.val < 0.001)',tmp$estimate)
          else
            main <- sprintf('Corr: %.3f (p.val = %.3f)',tmp$estimate, tmp$p.value)
        } 
        # plot
        par(mfrow=c(1,1), mar=c(4.1,4.1,1,1))
        xlim <- c(min(x)/(1+padx), max(x)*(1+padx))
        ylim <- c(min(y)/(1+pady), max(y)*(1+pady))
        plot(x, y, xlab=varx, ylab=vary, main=main, col=rgb(.3,.3,.3,.6), pch=16, xlim=xlim, ylim=ylim) 
        abline(v=mean(x),lty=2,col=rgb(.1,.1,.1,.4)); abline(h=mean(y),lty=2,col=rgb(.1,.1,.1,.4))
        contour(kde, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE,lwd=2)
      }
    },
    
    setEst = function(mcmc.output=NA, probs=NA, burninProportion=.2) {
      mcmc.output <- ifelse(is.na(mcmc.output), self$mcmc, mcmc.output)
      if(all(is.na(mcmc.output))) stop('mcmc.output must be already set or provided as argument')
      for (param in self$params) {
        samps <- unlist(sapply(mcmc.output,function(x){
          len <- nrow(x)
          burn <- ceiling(burninProportion*len)
          x[burn:len, param]
        }))
        self$est[[param]] <- quantile(samps, self$probs)
      }
    },
    
    getEst = function() {
      return(self$est)
    }
  )
)






