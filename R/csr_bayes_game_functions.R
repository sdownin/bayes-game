library(rjags)
library(ggplot2)
library(coda)
library(mcmcplots)
library(lattice)
library(latticeExtra)
library(plyr)

setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')

#--------------------------- FUNCTIONS -----------------------------------------

##
# 
# @argument x [list]; x$z can be n-length vector of {1,0}
# @returns [vector]  demand share for each x$z
##
share.list <- function(x,k=1) 
{
  out <- sapply(x$z, function(z){
    th1 <- x$p2 * (x$v1 + x$omega * x$sig1 * z)
    th2 <- x$p1 * (x$v2 + x$omega * x$sig2 * z)
    if(k==1) {
      num <- (th1/th2) * x$J1^x$rho
      denom <- (th1/th2) * x$J1^x$rho  + x$J2^x$rho
    } else {
      num <- (th2/th1) * x$J2^x$rho
      denom <- (th2/th1) * x$J2^x$rho  + x$J1^x$rho   
    }   
    return(num / denom)
  })
  return(out)
}

##
#
#
##
share <- function(p1,p2,v1,v2,sig1,sig2,J1,J2,omega,z,rho,k=1) 
{
  out <- sapply(z, function(z_i){
    th1 <- p2 * (v1 + omega * sig1 * z_i)
    th2 <- p1 * (v2 + omega * sig2 * z_i)
    if(k==1) {
      num <- (th1/th2) * J1^rho
      denom <- (th1/th2) * J1^rho  + J2^rho
    } else {
      num <- (th2/th1) * J2^rho
      denom <- (th2/th1) * J2^rho  + J1^rho   
    }   
    return(num / denom)
  })
  return(out)
}

##
# 
# @returns list to be used as `x` argument in demand share() function
##
getTheta <- function(v1 = 1, v2 = 1, J1 = 40, J2 = 60, p1 = 10, p2 = 10, 
                     sig1 = 1, sig2 = 0, omega = 1.5 , rho = 1, z = NA)
{
  return(list(v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,sig1=sig1,sig2=sig2,omega=omega,rho=rho,z=z))
}

##
#
#
##
getMaxPsi <- function(omega,rho,x,p1,p2,J1,J2)
{
  r1 <- x$r1; w1 <- x$w1; v1 <- x$v1; v2 <- x$v2;
  a <- (r1*p1-w1)*(omega/(v1+omega))
  b.numer <- (v2 + omega) * p1 * J2^rho
  b.denom <- v1*p2*J1^rho + (v2+omega)*p1*J2^rho
  return (a * b.numer / b.denom)
}

##
#
# @returns [data.frame] simulated data for evaluating CSR Bayes game via Gibbs sampling
##
getSimulation <- function(q=.5,Y=2000,N=1000, setSeed=TRUE,theta=NA)
{
  if(any(is.na(theta)))
    theta <- getTheta()
  L <- floor(Y/mean(theta$p1,theta$p2))  ## number of purchases each period per person i=1,...,N
  #
  if(setSeed)
    set.seed(111)
  df.sim <- data.frame( z=rbinom(n=N, size = 1, q))
  df.sim$s <- share.list(getTheta(z=df.sim$z))
  df.sim$G1 <- sapply(df.sim$s, function(x)rbinom(n=1, size = L, x))
  df.sim$G2 <- L - df.sim$G1 
  df.sim$sig1 <- theta$sig1
  df.sim$sig2 <- theta$sig2
  return(df.sim)
}

##
#
#
##
getCumuAvg <- function(mcmc.output, n.chains) 
{
  output <- mcmc.output
  df.cumavg <- data.frame(qavg=NA,x=NA,chain=NA)
  for (run in 1:n.chains) {
    x <- output[[run]]
    df.tmp <- data.frame(qavg=cumsum(x)/seq_len(length(x)), 
                         x=seq_len(length(x)), 
                         chain=paste0('chain ',run))
    df.cumavg <- na.omit( rbind(df.cumavg, df.tmp) )
  }
  return(df.cumavg)
}

##
# check if Markov Chain found stationary distribution with Autocorrelation  NHST
# @returns [data.frame]
##
autocorrTestMcmc <- function(mcmc.output, nlags=20, pvalOnly=FALSE, type='Ljung-Box')
{
  if (is.na(type) | type != 'Ljung-Box')
    type <- 'Box-Pierce'
  tests <- sapply(mcmc.output, function(x)Box.test(x, lag=nlags, type=type))
  out <- data.frame(t(tests))
  if(pvalOnly)
    return( c(p.value.avg=mean(sapply(out$p.value,function(x)x))) )
  return(out)
}

##
#
#
##
getGrowingVector <- function(N0,Tau,growth)
{
  out <- lapply(seq(0,Tau-1),function(t){
    n <- ceiling(N0*(1+growth)^t)
    return(rep(NA,n))
  })
  return(out)
}

##
#
#
##
getPsi <- function(gamma,y,p,B)
{
  return(gamma / ((y/p)*B))
}

##
#
#
##
getB <- function(s,m,b,d)
{
  newB <- s*m + b*(1-d)
  return(ifelse(newB > 0, newB, 0))
}

##
#
#
##
getJ <- function(y,rho,gamma,c,B,f,J,dj)
{
  netMargProfit <- ((rho+1)/rho) - (gamma/(rho*c))
  newJ <- y*netMargProfit*(B/f) + J*(1-dj)
  return(ifelse(newJ > 0, newJ, 0))
}

##
#
#
##
getG <- function(s,L,M,seed=1111)
{
  set.seed(seed)
  size <- min(M, length(s))
  s.sample <- sample(s,size,replace = F)
  return( sapply(s.sample, function(s)rbinom(n=1, size = L, s)) )
}

##
#
#
##
getQty <- function(y,p,netB,G)
{
  return(netB*(y/p) + G)
}

##
#
#
##
getPi <- function(r,omega,psi,rho,c,Q,O,y)
{
  netMargProfit <-  r - ( (omega+psi)/(rho*c) )
  return( netMargProfit*Q - (O/y) )
}

##
#
#
##
muBeta <- function(a,b)
{
  return(a/(a+b))
}

##
#
#
##
varBeta <- function(a,b)
{
  num <- a*b
  denom <- (a+b)^2 * (a+b+1)
  return(num/denom)
}

##
#
#
##
getModelstring <- function(sig1,sig2)
{
  paste0("model{
         for (i in 1:n) {
         z[i] ~ dbern(q)
         th1[i] <- p2*(v1 + omega*sig1*z[i])
         th2[i] <- p1*(v2 + omega*sig2*z[i])
         s[i] <- ",ifelse(sig2 > sig1,
                          "((th2[i]/th1[i])*pow(J2, rho)) / ( pow(J1, rho) + (th2[i]/th1[i])*pow(J2, rho) )",
                          "((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )"
         ),"
      G[i] ~ dbinom(s[i],L)
}
q ~ dbeta(h1t,h2t)
}")
}

##
#
#
##
getQhatMcmc <- function(data, modelstring, variables=c('q'), n.chains=2, n.adapt=3000,n.iter.update=3000, n.iter.samples=3000, thin=3, seed=1111)
{
  cat('Starting MCMC . . . \n')
  set.seed(1111)
  model <- jags.model(textConnection(modelstring),
                      data=data,n.adapt=n.adapt,n.chains=n.chains )
  ## BURN IN ITERATIONS TO APPROACH STATIONARY DISTRIBUTION
  update(model, n.iter=n.iter.update)
  ## ADAPTIVE SAMPLES FROM STATIONARY DISTRIBUTION TO SIMULATE POSTERIOR 
  output <- coda.samples(model=model, variable.names=variables,n.iter=n.iter.samples, thin=thin)
  return(output)
}

##
#
#
##
getQhatEst <- function(mcmc.output, probs, burninProportion=.2)
{
  samps <- unlist(sapply(mcmc.output,function(x){
    len <- length(x)
    burn <- ceiling(burninProportion*len)
    x[burn:len]
  }))
  return(quantile(samps, probs))
}

##
#
#
##
plotMCMCdiagnostics <- function(output,t,param='q')
{
  n.chains <- length(output)
  par(mfrow=c(2,2),mar=c(4,3,2,2))
  mcmcplots::denplot(output,style = 'plain',auto.layout = F,main=sprintf("Density of %s (t=%s)",param,t))
  mcmcplots::traplot(output,style = 'plain',auto.layout = F,main=sprintf("Trace of %s (t=%s)",param,t))
  mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main=sprintf("Thinned Running Mean of %s (t=%s)",param,t))
  mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main=sprintf("Autocorrelation of %s (t=%s)",param,t))
}

##
#
#
##
getNaiveMCMCseMean <- function(mcmc.output)
{
  su <- summary(mcmc.output)
  mu <- su$statistics['Mean']
  t.stat <- unname(mu/su$statistics['SD'])
  p.val <- pt(q = t.stat, lower.tail = F, df = su$end-su$start-1)
  return(print(sprintf('MCMC mean point estimate = %.4f (p-val = %.4f)',mu,p.val)))
}

##
#
#
##
getNIterRafDiag <- function(mcmc.samp.output, q=0.025, r=0.005, s=0.95)
{
  raf.test <- try(
    raftery.diag(c(mcmc.samp.output), q=q, r=r, s=s), 
    silent=T
  )
  if (inherits(raf.test, 'try-error'))
    stop(sprintf('the error level r=%s is too high. Try again with smaller r.',r))
  nmin <- as.numeric(raf.test$resmatrix[2])
  cat(sprintf('Raftery & Lewis N_min interations = %s',nmin))
  return(nmin)
}

##
#
#
##
getHeidConvDiag <- function(mcmc.output, eps=0.1, pvalue=0.05)
{
  cat('\nHeidelberg-Welch diagnostic:\n')
  h <- heidel.diag(mcmc.output, eps=eps, pvalue=pvalue)
  print(h)
  return(all(sapply(h, function(x)x[1,c('stest','htest')]==1)))
}

##
#
#
##
getQstarSig20 <- function(x,p1,p2,J1,J2,B1)
{
  omega <- x$omega; rho <- x$rho; r1 <- x$r1; c1 <- x$c1; w1 <- x$w1; v1 <- x$v1; v2 <- x$v2; 
  gamma1 <- x$gamma1; y <- x$Y
  Qty <- (y/p1) * B1
  psi1 <- gamma1 / Qty
  num <- psi1*v1*( p2*(v1+omega)*J1^rho  + p1*v2*J2^rho )
  denom <- omega*v2*p1*J2^rho * ( r1*rho*c1 - (w1+psi1) )
  return(num/denom)
}

##
#
#
##
getQstarSig21 <- function(x,p1,p2,J1,J2,B1)
{
  omega <- x$omega; rho <- x$rho; r1 <- x$r1; c1 <- x$c1; w1 <- x$w1; v1 <- x$v1; v2 <- x$v2; 
  gamma1 <- x$gamma1; y <- x$Y
  psi1 <- gamma1 / ((y/p1)*B1)
  num <- (psi1*v1*p2*J1^rho)/(v1*p2*J1^rho + v2*p1*J2^rho)
  den.a <- ((r1*rho*c1-(w1+psi1))*(v1+omega)*p2*J1^rho)/((v1+omega)*p2*J1^rho + v2*p1*J2^rho)
  den.b <- ((r1*rho*c1-w1)*v1*p2*J1^rho)/(v1*p2*J1^rho + (v2+omega)*p1*J2^rho)
  den.c <- (psi1*v1*p2*J1^rho)/(v1*p2*J1^rho + v2*p1*J2^rho)
  return(num/(den.a+den.b+den.c))
}

# ##
# #
# #
# ##
# getSigmaStar <- function(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1, qhat)
# {
#   ## FIND OPPONENT CUTOFF q
#   qstar2.if.sig1.0 <- getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
#   qstar2.if.sig1.1 <- getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
#   # FIND OPPONENT STRATEGY DECISION AS FUNTION OF OWN STRATEGY 
#   sig2.if.sig1.0 <- ifelse(qstar2.if.sig1.0 > qhat, 1, 0)
#   sig2.if.sig1.1 <- ifelse(qstar2.if.sig1.1 > qhat, 1, 0)
#   ## KNOWING WHAT THE OPPONENT WILL PLAY IF WE PLAY 0,  GET PROFITABLE STRATEGY COMPARED TO qhat
#   qstar1.of.sig2.if.sig1.0 <- ifelse(sig2.if.sig1.0==0,
#                                      getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1),
#                                      getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1))
#   sig1.of.sig2.if.sig1.0 <- ifelse(qstar1.of.sig2.if.sig1.0 > qhat, 1, 0)
#   ## KNOWING WHAT THE OPPONENT WILL PLAY IF WE PLAY 1,  GET PROFITABLE STRATEGY COMPARED TO qhat
#   qstar1.of.sig2.if.sig1.1 <- ifelse(sig2.if.sig1.1==0,
#                                      getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1),
#                                      getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1))
#   sig1.of.sig2.if.sig1.1 <- ifelse(qstar1.of.sig2.if.sig1.1 > qhat, 1, 0)
#   return()  # ??????????????????????????????????????
# }
getSigmaStar <- function()
{
  
}

##
# Main Function to Play CSR Bayes game
#   with MCMC Gibbs sampling to `learn` distribution of hedonic buyers (q)
# @param x [list] all game parameters
# @returns [list] game outcomes per variable and MCMC samples
##
playCsrBayesGame <- function(x, learn=TRUE)
{
  x$N <- ceiling(x$N0*(1+x$growth)^(x$Tau-1))
  ## allocate game array
  l <- list(
    q.true=x$q
    , M=rep(0,x$Tau)
    , qhat=list(mcmc=list(),
                est=data.frame(L99=rep(0,x$Tau),L95=rep(0,x$Tau),
                               mu=rep(0,x$Tau),
                               U95=rep(0,x$Tau),U99=rep(0,x$Tau)))
    , qstar=data.frame(qstar1=rep(0,x$Tau), qstar2=rep(0,x$Tau))
    , p=data.frame(p1=rep(10,x$Tau), p2=rep(10,x$Tau))
    , gamma=data.frame(gamma1=rep(0,x$Tau), gamma2=rep(0,x$Tau))
    , psi=data.frame(psi1=rep(0,x$Tau), psi2=rep(0,x$Tau))
    , f=data.frame(f1=rep(1,x$Tau), f2=rep(1,x$Tau))
    , O=data.frame(O1=rep(1,x$Tau), O2=rep(1,x$Tau))
    , J=data.frame(J1=rep(0,x$Tau),J2=rep(0,x$Tau))
    , B=data.frame(B1=rep(0,x$Tau),B2=rep(0,x$Tau))
    , h=data.frame(h1=rep(0,x$Tau),h2=rep(0,x$Tau))
    , sig=data.frame(sig1=rep(0,x$Tau), sig2=rep(0,x$Tau))
    , Q=data.frame(Q1=rep(0,x$Tau),Q2=rep(0,x$Tau))
    , Pi=data.frame(Pi1=rep(0,x$Tau),Pi2=rep(0,x$Tau))
    , z=list()  ##getGrowingVector(x$N0,x$Tau,x$growth)
    , s=list()  ##getGrowingVector(x$N0,x$Tau,x$growth)
    , G=list(G1=list(),G2=list())
  )
  
  t <- 1
  
  #--------------------------------- INITIAL PERIOD ------------------------------------------
  ## TEMP VALUES TO BE REPLACED BY LEARNING AFTER INITAL PERIOD
  l$qhat$est[t, ] <- quantile(rbeta(1e4, x$a1, x$a2), probs = x$probs)
  l$qhat$mcmc[t] <- NA
  l$qstar$qstar1[1] <- .5 ## ??????????????????
  l$qstar$qstar2[1] <- .5 ## ??????????????????
  
  ## Initial values
  l$J$J1[t] <- 50
  l$J$J2[t] <- 250
  l$B$B1[t] <- 500
  l$B$B2[t] <- 2500
  l$p$p1[t] <- x$c1 * x$rho
  l$p$p2[t] <- x$c2 * x$rho
  l$sig$sig1[t] <- 0
  l$sig$sig2[t] <- 0
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1,x$gamma2, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, getPsi(l$gamma$gamma2[t],x$Y,l$p$p2[t],l$B$B2[t]), 0)
  l$M[t] <- round(x$db1*l$B$B1[t] + x$db2*l$B$B2[t])
  l$z[[t]] <- rbinom(l$M[t], 1, x$q)  ##  CHANGING TO GROUND TRUE  q
  l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
  l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])
  
  
  
  # LIST demand share
  l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$omega, l$z[[t]],
                    x$rho, k=1)
  l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
  l$G$G2[[t]] <- getG( 1-l$s[[t]], l$L$L2[t], l$M[t])
  
  ## Qstar THRESHOLD
  # qstar1.0 <- getQstarSig20(x$omega,x$rho,x$r1,x$c1,x$w1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t])
  # qstar1.1 <- getQstarSig21(x$omega,x$rho,x$r1,x$c1,x$w1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t])
  # qstar2.0 <- getQstarSig20(x$omega,x$rho,x$r2,x$c2,x$w2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t])
  # qstar2.1 <- getQstarSig21(x$omega,x$rho,x$r2,x$c2,x$w2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t])
  # l$qstar$qstar1[t] <- getSigmaStar(x$omega,x$rho,x$r1,x$c1,x$w1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t], l$qhat$est$mu[t]) ## ??????????????????
  # l$qstar$qstar2[t] <- getSigmaStar(x$omega,x$rho,x$r2,x$c2,x$w2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t], l$qhat$est$mu[t]) ## ??????????????????
  
  ## LEARN Qhat
  if(learn) {
    data <- list(G=l$G$G1[t][[1]],  L=l$L$L1[t], 
                 n=l$M[t]-1, ## WHY NEED TO DROP ONE OBS??             ## NOTE NO Z HERE (UNOBSERVED)
                 sig1=l$sig$sig1[t], sig2=l$sig$sig2[t],
                 J1=l$J$J1[t],J2=l$J$J2[t],p1=l$p$p1[t],p2=l$p$p2[t],
                 v1=x$v1,v2=x$v2,
                 omega=ifelse(l$qhat$est$mu[t]>0, x$omega, 0),    ## ensure no signal when q=0
                 rho=x$rho,
                 h1t=x$a1 + l$h$h1[t], h2t=x$a2 + l$h$h2[t])
    modelstring1 <- getModelstring(l$sig$sig1[t], l$sig$sig2[t])
    l$qhat$mcmc[[t]]  <- getQhatMcmc(data, modelstring1, variables=c('q'), 
                                     n.chains=3, n.adapt=100,n.iter.update=100, 
                                     n.iter.samples=100, thin=2, seed=1111)
    l$qhat$est[t, ] <- getQhatEst(l$qhat$mcmc[[t]], probs=x$probs, burninProportion = .2)
    
    ## CHECK N_MIN  ITERATIONS WITH RAFTERY & LEWIS DIAGNOSTIC 
    ## ERROR:  +/- 1%
    x$n.iter <- getNIterRafDiag(l$qhat$mcmc[[t]], q = .025, r=.01, s=0.95)    
  }
  
  
  ## TEST LEARNING VIA AUTOCORRELATION
  # l$h$h1[t] <- x$ep * ( sum( sapply(seq_len(t),function(ii)sum(l$z[[ii]])) ) )
  # l$h$h2[t] <- x$ep * ( sum(sapply(seq_len(t), function(ii)length(l$z[[ii]])))  - sum( sapply(seq_len(t),function(ii)sum(l$z[[ii]])) ) )
  # ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=T, type='Ljung-Box')
  # if ( ac > x$learningThreshold) {
  #   l$h$h1[t] <- x$ep *  sum(l$z[[t]]) 
  #   l$h$h2[t] <- x$ep * ( length(l$z[[t]]) -  sum(l$z[[t]]) )
  # } else {
  #   l$h$h1[t] <- 0
  #   l$h$h2[t] <- 0
  # }
  ## USE CHI_SQUARE STATISTIC FROM AUTOCORR TEST TO DOWNWEIGHT EVIDENCE LEARNED FROM THIS PERIOD MCMC
  if(x$downweight & learn) {
    ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=F, type='Ljung-Box')
    l$h$h1[t] <- x$a1 + sum(l$z[[t]])/log(mean(unlist(ac$statistic)))
    l$h$h2[t] <- x$a2 + ( length(l$z[[t]]) - sum(l$z[[t]]) ) /log(mean(unlist(ac$statistic)))
  } else {
    l$h$h1[t] <- x$a1 + sum(l$z[[t]])
    l$h$h2[t] <- x$a2 + ( length(l$z[[t]]) - sum(l$z[[t]]) )
  }
  
  
  #### OUTCOME2
  ## Quantity
  l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
  l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
  ## Platform Operator Profit
  l$Pi$Pi1[t] <- getPi(x$r1,x$w1,l$psi$psi1[t],x$rho,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
  l$Pi$Pi2[t] <- getPi(x$r2,x$w2,l$psi$psi2[t],x$rho,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
  
  
  
  
  #---------------------------------- MAIN GAME LOOP ----------------------------------------------
  
  for (t in 2:x$Tau)
  {
    cat(paste0('\nt: ',t,'\n'))
    
    ## STRATEGY DECISION VARIABLES
    # l$qstar$qstar1[t] <- getSigmaStar(x$omega,x$rho,x$r1,x$c1,x$w1,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma1,l$B$B1[t-1], l$qhat$est$mu[t-1]) ## ??????????????????
    # l$qstar$qstar2[t] <- getSigmaStar(x$omega,x$rho,x$r2,x$c2,x$w2,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma2,l$B$B2[t-1], l$qhat$est$mu[t-1]) ## ??????????????????
    l$qstar$qstar1[t] <-  ifelse(l$sig$sig2[t]==0,
                                 getQstarSig20(x,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],l$B$B1[t-1]),
                                 getQstarSig21(x,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],l$B$B1[t-1])
    )
    l$qstar$qstar2[t] <- ifelse(l$sig$sig2[t]==1,
                                getQstarSig20(x,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],l$B$B2[t-1]),
                                getQstarSig21(x,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],l$B$B2[t-1])
    )
    ## CSR STRATEGIES
    if(!is.na(x$sig1.fixed)) { ## FIXED STRATEGY 
      l$sig$sig1[t] <- ifelse(!is.na(x$sig1.fixed[t]), x$sig1.fixed[t],  l$sig$sig1[t-1])
    } else {                     ## BAYES LEARNED STRATEGY
      l$sig$sig1[t] <- ifelse(t<4,0,1) ##ifelse(l$qhat$est$mu[t-1] > l$qstar$qstar1[t], 1, 0)
    }
    
    if(!is.na(x$sig2.fixed)) { ## FIXED STRATEGY 
      l$sig$sig2[t] <- ifelse(!is.na(x$sig2.fixed[t]), x$sig2.fixed[t],  l$sig$sig2[t-1])
    } else {                     ## BAYES LEARNED STRATEGY
      l$sig$sig2[t] <- ifelse(t<4,0,1) ##ifelse(l$qhat$est$mu[t-1] > l$qstar$qstar1[t], 1, 0)
    }
    
    ## UPDATE STRATEGY CHANGES
    if(l$sig$sig1[t] != l$sig$sig1[t-1] )
      x$t1.change <- c(x$t1.change, t)
    if(l$sig$sig2[t] != l$sig$sig2[t-1] )
      x$t2.change <- c(x$t2.change, t)
    
    ## CSR CONTINGENT COSTS
    l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1, x$gamma1, 0)
    l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
    l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, x$psi1, 0)
    l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, x$psi2, 0)
    
    ## PERIOD PARAMETERS
    l$p$p1[t] <- x$c1 * x$rho
    l$p$p2[t] <- x$c2 * x$rho
    l$f$f1[t] <- 1
    l$f$f2[t] <- 1
    l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
    l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])
    ## PREVIOS PERIOD DEPENDENT UPDATES
    m.temp <- round( x$db1*l$B$B1[t-1] + x$db2*l$B$B2[t-1] )
    l$M[t] <- ifelse( m.temp <= 1, 1, m.temp)
    l$J$J1[t] <- getJ(x$Y,x$rho,l$gamma$gamma1[t],x$c1,l$B$B1[t-1],l$f$f1[t],l$J$J1[t-1],x$dj1) 
    l$J$J2[t] <- getJ(x$Y,x$rho,l$gamma$gamma2[t],x$c2,l$B$B2[t-1],l$f$f2[t],l$J$J2[t-1],x$dj2)
    s1.avg <- sum(l$G$G1[[t-1]]) / (length(l$G$G1[[t-1]])*l$L$L1[t])
    s2.avg <- sum(l$G$G2[[t-1]]) / (length(l$G$G2[[t-1]])*l$L$L2[t])
    l$B$B1[t] <- getB(s1.avg, l$M[[t]], l$B$B1[t-1], x$db1)         
    l$B$B2[t] <- getB(s2.avg, l$M[[t]], l$B$B2[t-1], x$db2) 
    
    # SAMPLE MARKET ATTITUDES
    l$z[[t]] <- rbinom(l$M[t],1, x$q)
    
    # LIST demand share
    l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                      x$v1, x$v2,
                      l$sig$sig1[t], l$sig$sig2[t],
                      l$J$J1[t], l$J$J2[t],
                      x$omega, l$z[[t]],
                      x$rho, k=1)
    l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
    l$G$G2[[t]] <- getG(1-l$s[[t]], l$L$L2[t], l$M[t])
    
    
    ## LEARN Qhat
    if(learn) {
      data <- list(G=l$G$G1[t][[1]],  L=l$L$L1[t], 
                   n=l$M[t]-1,  ## NO Z HERE
                   sig1=l$sig$sig1[t], sig2=l$sig$sig2[t],
                   J1=l$J$J1[t],J2=l$J$J2[t],p1=l$p$p1[t],p2=l$p$p2[t],
                   v1=x$v1,v2=x$v2,
                   omega=ifelse(l$qhat$est$mu[t]>0, x$omega, 0),    ## ensure no signal when q=0
                   rho=x$rho,
                   h1t=x$a1 + l$h$h1[t-1], h2t=x$a2 + l$h$h2[t-1])
      modelstring1 <- getModelstring(l$sig$sig1[t], l$sig$sig2[t])
      l$qhat$mcmc[[t]] <- getQhatMcmc(data, modelstring1, variable=c('q'), 
                                      n.chains=3, n.adapt=1000,n.iter.update=x$n.iter, 
                                      n.iter.samples=x$n.iter, thin=3, seed=1111)
      l$qhat$est[t, ] <- getQhatEst(l$qhat$mcmc[[t]], probs=x$probs, burninProportion = .2)
      
      ## MCMC SAMPLE MEAN
      getNaiveMCMCseMean(l$qhat$mcmc[[t]])    
      
      ## MCMC DIAGNOSTICS
      getHeidConvDiag(l$qhat$mcmc[[t]])
      gelman.plot(l$qhat$mcmc[[t]], main=sprintf('\n\niteration: %s',t))
      plotMCMCdiagnostics(l$qhat$mcmc[[t]],t)      
    }
    
    
    ## TEST LEARNING VIA AUTOCORRELATION
    ## IF STRATEGIES ARE  THE SAME, AUTOCORRELATION SHOULD BE SIGNIFICANT --> DON'T COUNT THIS PERIOD z SAMPLE
    # ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=T, type='Ljung-Box')
    # if (  ac > x$learningThreshold ) {
    #   l$h$h1[t] <- x$ep*sum(l$z[[t]])  + l$h$h1[t-1]
    #   l$h$h2[t] <- x$ep*( length(l$z[[t]]) - sum(l$z[[t]]) ) + l$h$h2[t-1]
    # } else {
    #   l$h$h1[t] <- l$h$h1[t-1]
    #   l$h$h2[t] <- l$h$h2[t-1]
    # }
    ###
    ## USE CHI_SQUARE STATISTIC FROM AUTOCORR TEST TO DOWNWEIGHT EVIDENCE LEARNED FROM THIS PERIOD MCMC
    if(x$downweight & learn ) {
      ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=F, type='Ljung-Box')
      l$h$h1[t] <- l$h$h1[t-1] + sum(l$z[[t]])/log(mean(unlist(ac$statistic)))
      l$h$h2[t] <- l$h$h2[t-1] + ( length(l$z[[t]]) - sum(l$z[[t]]) ) /log(mean(unlist(ac$statistic)))
    } else {
      l$h$h1[t] <- l$h$h1[t-1] + sum(l$z[[t]])
      l$h$h2[t] <- l$h$h2[t-1] + ( length(l$z[[t]]) - sum(l$z[[t]]) ) 
    }
    
    #### OUTCOME2
    ## Quantity
    l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
    l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
    ## Platform Operator Profit
    l$Pi$Pi1[t] <- getPi(x$r1,x$w1,l$psi$psi1[t],x$rho,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
    l$Pi$Pi2[t] <- getPi(x$r2,x$w2,l$psi$psi2[t],x$rho,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
  }
  #--------------------------------- END MAIN GAME LOOP ------------------------------------------------
  return(l)
}

##
#
#
##
getCsrBayesGameSummaryPlots <- function(x,l)
{
  par(mfrow=c(3,2),mar=c(4,2.5,4,1))
  matplot(l$J/rowSums(l$J),xlab=('Game Period (t)'),ylim=c(0,1),type='o',pch=16,main=expression(tilde(J)));abline(v=c(x$t1.change,x$t2.change),lty=2)
  matplot(l$B/rowSums(l$B),xlab=('Game Period (t)'),ylim=c(0,1),type='o',pch=16,main=expression(tilde(B)));legend('topright',legend=1:2,col=1:2,lty=1:2,pch=16);abline(v=c(x$t1.change,x$t2.change),lty=2)
  matplot(l$Q/rowSums(l$Q),xlab=('Game Period (t)'),ylim=c(0,1),type='o',pch=16,main=expression(tilde(Qty)));abline(v=c(x$t1.change,x$t2.change),lty=2)
  matplot(l$Pi,            xlab=('Game Period (t)'),            type='o',pch=16,main=expression(pi));abline(h=0,col='black');abline(v=c(x$t1.change,x$t2.change),lty=2)
  matplot(l$h/rowSums(l$h),xlab=('Game Period (t)'),ylim=c(0,1),type='o',pch=16,main=expression(tilde(h)==Sigma[i]*z[i]),col=c('steelblue','darkgreen'))
  if(any(sapply(l$qhat$mcmc,function(x) length(na.omit(x)) > 0)))
    matplot(l$qhat$est,  xlab=('Game Period (t)'),ylim=c(0,1),type='o',pch=c(NA,NA,16,NA,NA),main=expression(hat(q)),lty=c(3,2,1,2,3),lwd=c(1,1,2,1,1),col=c('black','black','steelblue','black','black'));abline(h=x$q,col='black')
}

