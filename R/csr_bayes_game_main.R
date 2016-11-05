##
#
# CSR Bayes Game Main Simulation
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_functions_2.R'))
source(file.path(getwd(),'R','bvwg_r6_jags.R'))      # {.GLOBAL} JAGS



##
#
#
##
initCsrBayesGameConfig <- function(x)
{
  list( q.true=x$q
        , epsilon.true=x$epsilon
        , M=rep(0,x$Tau)
        , sim=list()
        , est=data.frame( L99=rep(0,x$Tau),L95=rep(0,x$Tau), mu=rep(0,x$Tau), U95=rep(0,x$Tau),U99=rep(0,x$Tau) )
        , p=data.frame(p1=rep(x$p1.0, x$Tau), p2=rep(x$p2.0, x$Tau))
        , gamma=data.frame(gamma1=rep(0,x$Tau), gamma2=rep(0,x$Tau))
        , psi=data.frame(psi1=rep(0,x$Tau), psi2=rep(0,x$Tau))
        , f=data.frame(f1=rep(1,x$Tau), f2=rep(1,x$Tau))
        , O=data.frame(O1=rep(1,x$Tau), O2=rep(1,x$Tau))
        , J=data.frame(J1=rep(0,x$Tau),J2=rep(0,x$Tau))
        , B=data.frame(B1=rep(0,x$Tau),B2=rep(0,x$Tau))
        , h=data.frame(h1=rep(0,x$Tau),h2=rep(0,x$Tau))
        , shape=rep(0,x$Tau)
        , rate=rep(0,x$Tau)
        , sig=data.frame(sig1=x$sig1, sig2=x$sig2)
        , Q=data.frame(Q1=rep(0,x$Tau),Q2=rep(0,x$Tau))
        , Pi=data.frame(Pi1=rep(0,x$Tau),Pi2=rep(0,x$Tau))
        , z=list()  
        , s=list() 
        , G=list(G1=list(),G2=list())
  )
}

##
# Main Function to Play CSR Bayes game
#   with MCMC Gibbs sampling to `learn` distribution of hedonic buyers (q)
# @param x [list] all game parameters
# @returns [list] game outcomes per variable and MCMC samples
##
playCsrBayesGame <- function(x, learn=TRUE, verbose=TRUE)
{
  l <- initCsrBayesGameConfig(x)
  t <- 1
  #--------------------------------- INITIAL PERIOD ------------------------------------------
  ## TEMP VALUES TO BE REPLACED BY LEARNING AFTER INITAL PERIOD
  l$est[t, ] <- quantile(rbeta(1e4, x$a1, x$a2), probs = x$probs)
  l$mcmc[[t]] <- NA
  
  ## Initial values
  l$J$J1[t] <- x$J1.0
  l$J$J2[t] <- x$J2.0
  l$B$B1[t] <- 4 * l$J$J1[t]
  l$B$B2[t] <- 4 * l$J$J2[t]
  l$p$p1[t] <- getPrice(1, x)
  l$p$p2[t] <- getPrice(2, x)
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1,x$gamma2, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, getPsi(l$gamma$gamma2[t],x$Y,l$p$p2[t],l$B$B2[t]), 0)
  l$M[t] <- round(x$db1*l$B$B1[t] + x$db2*l$B$B2[t])
  l$z[[t]] <- rbinom(l$M[t], 1, l$q.true)  ##  CHANGING TO GROUND TRUTH  q
  l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
  l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])
  
  # LIST demand share
  l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$omega, l$z[[t]],
                    x$epsilon, k=1)
  l$s.base <- share.base(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$omega, l$z[[t]],
                    x$epsilon, k=1)
  l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
  l$G$G2[[t]] <- getG( 1-l$s[[t]], l$L$L2[t], l$M[t])
  
  ## LEARN Qhat
  if(learn) {
    l <- learnBayesParams(x,l,t)
    ## CHECK N_MIN  ITERATIONS WITH RAFTERY & LEWIS DIAGNOSTIC 
    ## ERROR:  +/- 1%
    x$n.iter <- getNIterRafDiag(l$sim[[t]]$mcmc, q = .025, r=.01, s=0.95)    
  }
  
  #### OUTCOME2
  ## Quantity
  l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
  l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
  ## Platform Operator Profit
  l$Pi$Pi1[t] <- getPi(x$r1,x$w1,l$psi$psi1[t],x$epsilon,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
  l$Pi$Pi2[t] <- getPi(x$r2,x$w2,l$psi$psi2[t],x$epsilon,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
  
  #---------------------------------- MAIN GAME LOOP ------------------------------
  for (t in 2:x$Tau) {
    if(verbose) cat(paste0('\nt: ',t,'\n'))
    l <- updateGame(l, t)
  }
  #--------------------------------- END MAIN GAME LOOP ---------------------------
  return(l)
}

##
#
##
getModelStr <- function() {
  "model{
    for (i in 1:n) {
      z[i] ~ dbern(q)
      th1[i] <- p2*(v1 + omega*sig1*z[i])
      th2[i] <- p1*(v2 + omega*sig2*z[i])
      s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
      G[i] ~ dbinom( s[i], L )
    }
    epsilon ~ dgamma(shapet,ratet)
    q ~ dbeta(h1t,h2t)
  }"
}

##
#
#
##
learnBayesParams <- function(x,l,t)
{ 
  s.diff <- mean(l$s.base) - mean(l$s[[t]])
  s.cumu <- sum(unlist(l$G$G1)) / sum(unlist(c(l$G$G1,l$G$G2)))
  l$h$h1[t] = ifelse(t==1,x$a1,l$h$h1[t-1]) + ifelse(s.diff > 0, 1, 0)
  l$h$h2[t] = ifelse(t==1,x$a2,l$h$h2[t-1]) + ifelse(s.diff > 0, 1, 0)
  l$shape[t] = ifelse(t==1,1,l$shape[t-1])  + 1-s.cumu
  l$rate[t] = ifelse(t==1,1,l$rate[t-1])    + s.cumu
  ##
  data <- list(G=l$G$G1[t][[1]],  
               L=l$L$L1[t], 
               n=l$M[t]-1,  ## NO Z HERE
               sig1=l$sig$sig1[t], 
               sig2=l$sig$sig2[t],
               J1=l$J$J1[t],
               J2=l$J$J2[t],
               p1=l$p$p1[t],
               p2=l$p$p2[t],
               v1=x$v1,
               v2=x$v2,
               omega= ifelse(x$q==0,0,x$omega),    ## ensure no signal when q=0
               h1t=l$h$h1[t],
               h2t=l$h$h2[t],
               shapet= l$shape[t],
               ratet=l$rate[t] )
  ##
  l$sim[[t]] <- JAGS$new(getModelStr(), c('q','epsilon'))
  l$sim[[t]]$run(data, period=t)
  l$sim[[t]]$bivarPlot('q','epsilon', chain=1)
  l$est[t, ] <- getQhatEst(l$sim[[t]]$mcmc, probs=x$probs, burninProportion = .2)
  return(l)
}

##
#
#
##
updateGame <- function(l, t)
{
  ## CSR CONTINGENT COSTS
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1, x$gamma1, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, x$psi1, 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, x$psi2, 0)
  
  ## PERIOD PARAMETERS
  l$p$p1[t] <- getPrice(1, x)
  l$p$p2[t] <- getPrice(2, x)
  l$f$f1[t] <- 1
  l$f$f2[t] <- 1
  l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
  l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])
  ## PREVIOS PERIOD DEPENDENT UPDATES
  m.temp <- round( x$db1*l$B$B1[t-1] + x$db2*l$B$B2[t-1] )
  l$M[t] <- ifelse( m.temp <= 1, 1, m.temp)
  l$J$J1[t] <- getJ(x$Y,x$epsilon,l$gamma$gamma1[t],x$c1,l$B$B1[t-1],l$f$f1[t],l$J$J1[t-1],x$dj1) 
  l$J$J2[t] <- getJ(x$Y,x$epsilon,l$gamma$gamma2[t],x$c2,l$B$B2[t-1],l$f$f2[t],l$J$J2[t-1],x$dj2)
  s1.avg <- sum(l$G$G1[[t-1]]) / (length(l$G$G1[[t-1]])*l$L$L1[t])
  s2.avg <- sum(l$G$G2[[t-1]]) / (length(l$G$G2[[t-1]])*l$L$L2[t])
  l$B$B1[t] <- getB(s1.avg, l$M[[t]], l$B$B1[t-1], x$db1)         
  l$B$B2[t] <- getB(s2.avg, l$M[[t]], l$B$B2[t-1], x$db2) 
  
  # SAMPLE MARKET ATTITUDES
  l$z[[t]] <- rbinom(l$M[t],1, l$q.true)
  
  # LIST demand share
  l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$omega, l$z[[t]],
                    x$epsilon, k=1)
  l$s.base <- share.base(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$omega, l$z[[t]],
                    x$epsilon, k=1)
  l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
  l$G$G2[[t]] <- getG(1-l$s[[t]], l$L$L2[t], l$M[t])
  
  ## LEARN Qhat
  if(learn) {
    l <- learnBayesParams(x,l,t)
  }
  
  #### OUTCOME2
  ## Quantity
  l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
  l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
  ## Platform Operator Profit
  l$Pi$Pi1[t] <- getPi(x$r1,x$w1,l$psi$psi1[t],x$epsilon,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
  l$Pi$Pi2[t] <- getPi(x$r2,x$w2,l$psi$psi2[t],x$epsilon,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
  
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
}

