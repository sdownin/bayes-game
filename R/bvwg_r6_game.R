library(R6)
#
source(file.path(getwd(),'R','bvwg_r6_jags.R'))      # {.GLOBAL} JAGS
source(file.path(getwd(),'R','bvwg_r6_platform.R'))  # {.GLOBAL} PLATFORM


GAME <- R6Class("GAME", 
  public = list(
    name = NA,
    pl = list(),     # players list
    t = 1,      # game periods
    l = list(        # config list
        q.true = .3
      , epsilon.true = .9
      , omega=1.5
      , Y=1000
      , ep=1e-1
      , N0=50
      , Tau=10
      , probs=c(.005,.025,.5,.975,.995)
      , learningThreshold=.05
      , n.iter=1000
      , downweight=TRUE
    ),
    
    initialize = function(players=c('a','b'), modelString='', params=c(), name=NA, l=list()) {
      self$name <- name
      for (i in 1:length(players) ) {
        self$pl <- c(self$pl, PLATFORM$new(name=players[i]) )
      }
      self$set(l)
      self$setup(modelString, params)
    },
    
    add = function(l) {  self$l <- c(self$l, l) },
    set = function(l) {
      if (typeof(l) == 'list') {
        for (item in names(l)) {
          self$l[[item]] <- l[[item]]
        }
      }
    },
    
    setup = function(modelString, params) {
      t <- self$l$Tau
      self$l$M=rep(0,t)
      self$l$qhat=list(mcmc = JAGS$new(modelString, params),
                       est = data.frame(L99=rep(0,t),L95=rep(0,t),mu=rep(0,t),U95=rep(0,t),U99=rep(0,t)))
      self$l$z=list()  ##getGrowingVector(x$N0,x$Tau,x$growth)
      self$l$s=list()  ##getGrowingVector(x$N0,x$Tau,x$growth)
      t <- self$l$Tau
      for (i in 1:length(self$pl)) {
        self$pl[[i]]$l$p <- rep(10, t)
        self$pl[[i]]$l$gamma  <- rep(0, t)
        self$pl[[i]]$l$psi <- rep(0, t)
        self$pl[[i]]$l$f <- rep(1, t)
        self$pl[[i]]$l$O <- rep(1, t)
        self$pl[[i]]$l$J <- rep(0, t)
        self$pl[[i]]$l$B <- rep(0, t)
        self$pl[[i]]$l$h <- rep(0, t)
        self$pl[[i]]$l$sig <- rep(0, t)
        self$pl[[i]]$l$Q <- rep(0, t)
        self$pl[[i]]$l$Pi <- rep(0, t)
        self$pl[[i]]$l$G <- list()
        ## Initial values
        self$pl[[i]]$J[t] <- 20
        self$pl[[i]]$B[t] <- 8 * l$J$J1[t]
        self$pl[[i]]$p[t] <- getPrice(1, x)
        self$pl[[i]]$sig[t] <- 0
        self$pl[[i]]$gamma[t] <- ifelse(self$pl[[i]]$sig[t]==1,self$l$g, 0)
        l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
        l$M[t] <- round(x$db1*l$B$B1[t] + x$db2*l$B$B2[t])
        l$z[[t]] <- rbinom(l$M[t], 1, l$q.epsilon[t])  ##  CHANGING TO GROUND TRUTH  q
        l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
      }
      self$l$q.noisy <- sapply(rnorm(self$l$Tau, mean = self$l$q.true, sd = .05), function(x){
        ifelse(x < 0, 0, ifelse(x > 1, 1, x))
      })
    },
    
    increment = function() {
      
    }
    
  )
)












#--------------------------------- INITIAL PERIOD ------------------------------------------

## Initial values
l$J$J1[t] <- 20
l$B$B1[t] <- 8 * l$J$J1[t]
l$p$p1[t] <- getPrice(1, x)
l$sig$sig1[t] <- 0
l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1,x$gamma2, 0)
l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
l$M[t] <- round(x$db1*l$B$B1[t] + x$db2*l$B$B2[t])
l$z[[t]] <- rbinom(l$M[t], 1, l$q.epsilon[t])  ##  CHANGING TO GROUND TRUTH  q
l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])



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



#------------------------------------------------------------------------

g1 <- GAME$new()


## MAIN GAME CALL
## SET STRATEGY
x$t1.change <- 1
x$t2.change <- 4
x$sig1.fixed <- c(rep(0,x$t1.change),rep(1,x$Tau-x$t1.change))
x$sig2.fixed <- c(rep(0,x$t2.change),rep(1,x$Tau-x$t2.change))