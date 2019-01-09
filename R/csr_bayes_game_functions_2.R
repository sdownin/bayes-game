# library(rjags)
library(ggplot2)
library(coda)
# library(mcmcplots)
library(lattice)
library(latticeExtra)
library(plyr)
library(parallel)
# library(fields)
# library(random)
# library(foreach)
# library(snow)

setwd('C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')


#--------------------------- FUNCTIONS -----------------------------------------
##
#
##
ciPlot <- function(x, ...) {
  len <- ncol(x)
  matplot(x, type='o', pch=c(NA,NA,16,NA,NA), col=c(1,1,2,1,1), lty=c(3,2,1,2,3), lwd=c(1,1,2,1,1), ...)
}

##
#
##
isGoodClustering <- function(cl, cutoff=0.85, verbose=TRUE)
{
  betweenSSportion <- cl$betweenss / (sum(cl$withinss) + cl$betweenss)
  if(verbose)
    cat(sprintf('\nbetween SumSq portion: %.3f\n', betweenSSportion))
  return( cutoff < betweenSSportion )
}

##
#
##
getBetweenSS <- function(cl)
{
  return(cl$betweenss / (sum(cl$withinss) + cl$betweenss))
}

##
#
##
getDemandClusters <- function(l, t, k=2, maxcl=5, force=TRUE, all=FALSE)
{
  x <- l$G$G1[[t]]
  if( !force & l$sig$sig1[t]==l$sig$sig2[t]) {
    cat('\nSame CSR strategies. returning 1 cluster\n')
    return( kmeans(x = x, centers = 1)  )
  }
  is <- 1:maxcl
  cl <- list()
  if (typeof(k)=='character' | is.na(k)) {
    for (i in is) {
      cl[[i]] <- kmeans(x = x, centers = i)     
    }
    wss <- sapply(cl,function(x)sum(x$withinss))
    k <- which.min(diff(wss)) + 1
    if (all)
      return(cl)
    else
      return(cl[[k]])
  }
  return(kmeans(x = x, centers = k))
}

##
#
##
## see demand groups
plotDemandGroups <- function(x,cl) {
  hist(x, col='lightgray', breaks=15,main='Demand Groups'); abline(v=cl$centers, col='red', lwd=2)
}
##
# Get customer type vector Z = z_i in {0,1}
# @param vec          x    The observed demand 
# @param [[kmeans]]   cl   An object of class kmeans, kmeans clustering output
## 
getZ <- function(x,cl,t) {
  if (x$sig1[t] < x$sig2[t]) {
    z1index <- which.min(cl$centers)  
    z0index <- which.max(cl$centers)
  } else if (x$sig1[t] > x$sig2[t])  {
    z1index <- which.max(cl$centers)  
    z0index <- which.min(cl$centers)
  }
  Z <- sapply(cl$cluster, function(x)ifelse(x==z1index,1,0))  
}


##
#
##
getQfromEpsilonS <- function(epsilon,s,sig1=0,sig2=1,omega=1,J1=2,J2=20,v1=1,v2=1,p1=10,p2=10)
{
  s1 <- s
  s2 <- 1-s
  x1 <- (v1 * J1^epsilon) * p2 * s2
  x2 <- (v2 * J2^epsilon) * p1 * s1
  y1 <- (sig1 * J1^epsilon) * p2 * s2
  y2 <- (sig2 * J2^epsilon) * p1 * s1
  z <- (x1 - x2) / (y2 - y1)
  return((1/omega) * z)
}
##
#
##
persp.withcol <- function(x,y,z,pal,nb.col,...,xlg=TRUE,ylg=TRUE)
{
  colnames(z) <- y
  rownames(z) <- x
  
  nrz <- nrow(z)
  ncz <- ncol(z) 
  
  color <- pal(nb.col)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nb.col)
  par(xlog=xlg,ylog=ylg)
  
  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))
  s <- as.matrix(z)
  
  per <- persp(
    as.numeric(rownames(z)),
    as.numeric(colnames(z)),
    s,
    xlab=xname,ylab=yname,
    col=color[facetcol],
    ...
  )
}
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
      num <- (th1/th2) * x$J1^x$epsilon
      denom <- (th1/th2) * x$J1^x$epsilon  + x$J2^x$epsilon
    } else {
      num <- (th2/th1) * x$J2^x$epsilon
      denom <- (th2/th1) * x$J2^x$epsilon  + x$J1^x$epsilon   
    }   
    return(num / denom)
  })
  return(out)
}

##
#
#
##
share <- function(p1,p2,v1,v2,sig1,sig2,J1,J2,omega,z,epsilon,k=1) 
{
  out <- sapply(z, function(z_i){
    th1 <- p2 * (v1 + omega * sig1 * z_i)
    th2 <- p1 * (v2 + omega * sig2 * z_i)
    if(k==1) {
      num <- (th1/th2) * J1^epsilon
      denom <- (th1/th2) * J1^epsilon  + J2^epsilon
    } else {
      num <- (th2/th1) * J2^epsilon
      denom <- (th2/th1) * J2^epsilon  + J1^epsilon   
    }   
    return(num / denom)
  })
  return(out)
}

share.base <- function(p1,p2,v1,v2,sig1,sig2,J1,J2,omega,z,epsilon,k=1) 
{
  out <- sapply(z, function(z_i){
    th1 <- p2 * (v1)
    th2 <- p1 * (v2)
    if(k==1) {
      num <- (th1/th2) * J1^epsilon
      denom <- (th1/th2) * J1^epsilon  + J2^epsilon
    } else {
      num <- (th2/th1) * J2^epsilon
      denom <- (th2/th1) * J2^epsilon  + J1^epsilon   
    }   
    return(num / denom)
  })
  return(out)
}


##
#
##
expectedShare <- function(p1,p2,v1,v2,gamma1,gamma2,J1,J2,omega,epsilon,zH,zU,N,k=1) {
  sH <- .share(p1,p2,v1,v2,gamma1,gamma2,J1,J2,omega,epsilon,'H',k)
  sU <- .share(p1,p2,v1,v2,gamma1,gamma2,J1,J2,omega,epsilon,'U',k)
  return( sH*(zH / N) + sU*(zU / N) )
}
.share <- function(p1,p2,v1,v2,gamma1,gamma2,J1,J2,omega,epsilon,z,k=1) {
  gammak <- ifelse(k==1, gamma1, gamma2)
  pk <- ifelse(k==1, p2, p1)  ## opposite
  vk <- ifelse(k==1, v1, v2)
  Jk <- ifelse(k==1, J1, J2)
  ##
  zval <- ifelse( !is.character(z), z, ifelse(z=='H',1,0) )
  ## thetas
  th1 <- p2 * (vk + omega*gamma1*zval)
  th2 <- p1 * (vk + omega*gamma2*zval)
  thk <- pk * (vk + omega*gammak*zval)
  return( (thk * Jk^epsilon) / ((th1 * J1^epsilon)+(th2 * J2^epsilon)) )
}


##
# 
# @returns list to be used as `x` argument in demand share() function
##
getTheta <- function(v1 = 1, v2 = 1, J1 = 40, J2 = 60, p1 = 10, p2 = 10, 
                     sig1 = 1, sig2 = 0, omega = 1.5 , epsilon = 1, z = NA)
{
  return(list(v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,sig1=sig1,sig2=sig2,omega=omega,epsilon=epsilon,z=z))
}

##
#
#
##
getPlatformCsrCostPoolSepEq <- function(omega,epsilon,x,p1,p2,J1,J2)
{
  r1 <- x$r1; w1 <- x$w1; v1 <- x$v1; v2 <- x$v2;
  a <- (r1*p1-w1)*(omega/(v1+omega))
  b.numer <- (v2 + omega) * p1 * J2^epsilon
  b.denom <- v1*p2*J1^epsilon + (v2+omega)*p1*J2^epsilon
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
getJ <- function(y,epsilon,gamma,c,B,f,J,dj)
{
  net.eff <- epsilon + 1
  netMargProfit <- ((net.eff + 1)/net.eff) - (gamma/(net.eff*c))    
  newJ <- y*netMargProfit*(B/f) + J*(1-dj)
  return(ifelse(newJ > 0, newJ, 0))
}

##
#
#
##
getG <- function(s,L,M,seed=1111, dist='binomial')
{
  set.seed(seed)
  size <- min(M, length(s))
  s.sample <- sample(s,size,replace = F)
  G <- NA
  if (dist=='binomial')
    G <-  sapply(s.sample, function(s)rbinom(n=1, size = L, s))
  else if (dist == 'poisson')
    G <- sapply(s.sample, function(s)rpois(1, s*L))
  return(G)
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
getPi <- function(r,omega,psi,epsilon,c,Q,O,y)
{
  netMargProfit <-  r - ( (omega+psi)/(epsilon*c) )
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

# ##
# #
# #
# ##
# getModelstring <- function(sig1,sig2)
# {
#   paste0("model{
#          for (i in 1:n) {
#          z[i] ~ dbern(q)
#          th1[i] <- p2*(v1 + omega*sig1*z[i])
#          th2[i] <- p1*(v2 + omega*sig2*z[i])
#          s[i] <- ",ifelse(sig2 > sig1,
#                           "((th2[i]/th1[i])*pow(J2, rho)) / ( pow(J1, rho) + (th2[i]/th1[i])*pow(J2, rho) )",
#                           "((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )"
#          ),"
#       G[i] ~ dbinom(s[i],L)
# }
# q ~ dbeta(h1t,h2t)
# }")
# }

# ##
# #
# #
# ##
# getModelstring <- function(sig1,sig2)
# {
#   sprintf("
#     model{
#        for (i in 1:n) {
#           z[i] ~ dbern(q)
#           th1[i] <- p2*(v1 + omega*sig1*z[i])
#           th2[i] <- p1*(v2 + omega*sig2*z[i])
#           s[i] <- ((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )
#           G[i] ~ dbinom(s[i],L)
#         }
#       q ~ dbeta(h1t,h2t)
#     }
#   ")
# }

# getModelStr <- function() {
#   "model{
#   for (i in 1:n) {
#   z[i] ~ dbern(q)
#   th1[i] <- p2*(v1 + omega*sig1*z[i])
#   th2[i] <- p1*(v2 + omega*sig2*z[i])
#   s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
#   G[i] ~ dbinom( s[i], L )
#   }
#   epsilon ~ dgamma(shapet,ratet)
#   q ~ dbeta(h1t,h2t)
# }"
# }

# ##
# #
# #
# ##
# getModelstring2 <- function()
# {
#   sprintf("
#           model{
#           for (i in 1:n) {
#           z[i] ~ dbern(q)
#           th1[i] <- p2*(v1 + omega*sig1*z[i])
#           th2[i] <- p1*(v2 + omega*sig2*z[i])
#           s[i] <- ((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )
#           G[i] ~ dbinom(s[i],L)
#           }
#           q ~ dbeta(h1t,h2t)
#           rho ~ dnorm(mu,tau)
#           mu
#           }
#           ")
# }

# getModelStr <- function() {
#   "model{
#   for (i in 1:n) {
#   z[i] ~ dbern(q)
#   th1[i] <- p2*(v1 + omega*sig1*z[i])
#   th2[i] <- p1*(v2 + omega*sig2*z[i])
#   X[i] <- dbern( ((th1[i]/th2[i])*pow(J1, epsilon)) / ((th1[i]/th2[i])*pow(J1, epsilon)+pow(J2, epsilon)) )
#   }
#   epsilon ~ dgamma(shapet,ratet)
#   q ~ dbeta(h1t,h2t)
# }"
# }

##
# Gibbs Sampling estimate of q_hat
# JAGS implementation
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
# 
# n.iter <- 10000
# model.file <- "csr_bayes_game.bug"
# capture.output(cat(stringr::str_replace_all(getModelStr(),"\n","")), file=model.file)
# jags.inits <- function() {
#   return(list("q"=runif(1)))
# }
# output.mcmc <- as.mcmc(do.call(jags.parallel,
#                list(data = data,  parameters.to.save = c('q'),
#                     n.chains=4, n.iter=n.iter, n.burnin=ceiling(n.iter * .4),
#                     model.file=model.file)))
# out.mcmc <- as.mcmc(out)
# densplot(out.mcmc)
# 
# 
# # runJags <- function(data, modelstring, variables=c('q'), n.chains=2, n.adapt=3000,n.iter.update=3000, n.iter.samples=3000, thin=3, seed=1111)
# # {
#   cat('Starting MCMC . . . \n')
#   set.seed(1111)
#   clust <- makeCluster(4, type = "SOCK")
#   model <- jags.model(textConnection(getModelStr()),
#                       data=data,n.chains=4 )
#   ## BURN IN ITERATIONS TO APPROACH STATIONARY DISTRIBUTION
#   update(model, n.iter=1000)
#   ## ADAPTIVE SAMPLES FROM STATIONARY DISTRIBUTION TO SIMULATE POSTERIOR 
#   output <- coda.samples(model=model, variable.names='q',n.iter=1000, thin=2)
# #   return(output)
# # }

##
#
#
##
getParamEst <- function(mcmc.output, param, probs, burninProportion=.2)
{
  samps <- unlist(sapply(mcmc.output,function(x){
    len <- nrow(x)
    burn <- ceiling(burninProportion*len)
    x[burn:len, param]
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
  try(mcmcplots::denplot(output,style = 'plain',auto.layout = F,main=sprintf("Density of %s (t=%s)",param,t)))
  try(mcmcplots::traplot(output,style = 'plain',auto.layout = F,main=sprintf("Trace of %s (t=%s)",param,t)))
  try(mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main=sprintf("Thinned Running Mean of %s (t=%s)",param,t)))
  try(mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main=sprintf("Autocorrelation of %s (t=%s)",param,t)))
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
  cat(sprintf('Raftery & Lewis N_min interations = %s\n',nmin))
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
getPrice <- function(k, x, t)
{
  beta <- 1 + x$epsilon
  if(k %in% c(1,'1'))
    return( beta*(x$c1 + x$sig1[t]*x$gamma1 ) ) # + x$sig1[t]*x$phi1
  if(k %in% c(2,'2'))
    return( beta*(x$c2 + x$sig2[t]*x$gamma2 ) )
}

##
#
##
# getPriceFromCsrMarkup <- function(sig,epsilon,gamma,phiMarkup, c=1)
# {
#   basePrice <- (1+epsilon)*(c)
#   phi <- phiMarkup*basePrice / (1+epsilon)
#   return((1+epsilon)*(c + sig*gamma + sig*phi))
# }

getPriceFromCsrMarkup <- function(epsilon,gamma, c=1, verbose=FALSE)
{
  price <- (1+epsilon)*(c + gamma)
  if(verbose) cat(sprintf('price from CSR markup = %.3f\n',price))
  return(price)
}

getPriceFromCsrAndDiscountAction <- function(epsilon,gamma,discount, c=1, verbose=FALSE)
{
  discount.safe <- min(discount, epsilon)
  price <- (1+epsilon-discount.safe)*(c + gamma)
  if(verbose) cat(sprintf('PRICE(gamma,discount) = %.3f; discount > epsilon (warning?): %s\n',price, discount>epsilon))
  return(price)
}

# getPriceFromCsrMarkup <- function(sig,epsilon,gamma,phiMarkup, c=1, verbose=FALSE)
# {
#   ## IGNORE GAMMA USE PHIMARKUP AS COST `psi`
#   p0 <- c
#   basePrice <- p0 * (1 + phiMarkup)
#   price <- basePrice * (1+epsilon)
#   # if(verbose) cat(sprintf('phi = %.3f\n',(price-p0)/(1+epsilon)))
#   return(price)
# }

##
#
#
##
getMaxPsi <- function(k, x)
{
  epsilonTilde <- 1 + x$epsilon
  if(k %in% c(1,'1'))
    return(x$r1*epsilonTilde*x$c1 - x$w1)
  if(k %in% c(2,'2'))
    return(x$r2*epsilonTilde*x$c2 - x$w2)
}

