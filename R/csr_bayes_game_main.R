##
#
# CSR Bayes Game Main Simulation
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')
source(file.path(getwd(),'R','csr_bayes_game_functions_2.R'))
source(file.path(getwd(),'R','bvwg_r6_jags.R'))      # {.GLOBAL} JAGS

library(runjags)
library(random)
load.runjagsmodule()


##
# Main Function to Play CSR Bayes game
#   with MCMC Gibbs sampling to `learn` distribution of hedonic buyers (q)
# @param   [list] x  All game parameters
# @returns [list]    game outcomes per variable and MCMC samples
##
playCsrBayesGame <- function(x, ..., verbose=TRUE)
{
  l <- initCsrBayesGameConfig(x)
  t <- l$t
  #--------------------------------- INITIAL PERIOD ------------------------------------------
  ## LEARN Qhat
  if(x$learn) {
    l <- learnBayesParams(x,l,t)
    ## CHECK N_MIN  ITERATIONS WITH RAFTERY & LEWIS DIAGNOSTIC 
    ## ERROR:  +/- 1%
    r <- 0.005; counter <- 1
    l$n.iter <- try(getNIterRafDiag(l$sim[[t]]$mcmc, q = .025, r=r, s=0.95), silent = T)
    while ("try-error" %in% class(l$n.iter) & counter < 7) {
      cat(sprintf('Raftery & Lewis Diag redo: r = %.3f\n', r))
      r <- r*.5
      l$n.iter <- try(getNIterRafDiag(l$sim[[t]]$mcmc, q = .025, r=r, s=0.95), silent = T)
      counter <- counter + 1
    }
    l$n.iter <- ifelse("try-error" %in% class(x$n.iter), x$n.iter, l$n.iter)
    cat(sprintf('\nRaftery & Lewis Diag: %s\n', l$n.iter))
    ## RERUN with min samples 
    l <- learnBayesParams(x,l,t)
  }
  #### OUTCOME2
  ## Quantity
  l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
  l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
  ## Platform Operator Profit
  l$Pi$Pi1[t] <- getPi(x$r1,x$w1,l$psi$psi1[t],x$epsilon,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
  l$Pi$Pi2[t] <- getPi(x$r2,x$w2,l$psi$psi2[t],x$epsilon,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
  
  #---------------------------------- MAIN GAME LOOP ------------------------------
  if (x$Tau > 1) {
    for (t in 2:x$Tau) {
      l$t <- t
      if(verbose) cat(paste0('\nt: ',t,'\n'))
      l <- updateGame(x, l, t)
    }    
  }
  #--------------------------------- END MAIN GAME LOOP ---------------------------
  ## PLOTS
  if (x$learn) {
    for (param in x$params) {
      true <- switch(param, 'q'=l$q, 'epsilon'=l$epsilon)
      ciPlot(l$est[[param]], main=sprintf('MCMC Estimate for %s',param)); abline(v=x$t1.change);abline(h=true,col='gray')
    }
    ##cumuplot(na.omit(l$q.hat), main="q_hat Estimated from Inferred Z vector")
    hsummary <- sapply(1:nrow(l$h), function(i)quantile(rbeta(1000,l$h[i,1],l$h[i,2]),probs =l$probs))
    ciPlot(t(hsummary), main="q Beta Posterior from Inferred Buyer Types (Z)",xlab="period"); abline(v=c(x$t1.change,x$t2.change));abline(h=l$q,col='gray')
  }
  return(l)
}

##
#
#
##
initCsrBayesGameConfig <- function(x)
{
  l <- list( q=x$q
             , epsilon=x$epsilon
             , M=rep(0,x$Tau)
             , sim=list()
             , p=data.frame(p1=rep(NA, x$Tau), p2=rep(NA, x$Tau))
             , gamma=data.frame(gamma1=rep(0,x$Tau), gamma2=rep(0,x$Tau))
             , psi=data.frame(psi1=rep(0,x$Tau), psi2=rep(0,x$Tau))
             , f=data.frame(f1=rep(1,x$Tau), f2=rep(1,x$Tau))
             , O=data.frame(O1=rep(1,x$Tau), O2=rep(1,x$Tau))
             , J=data.frame(J1=rep(0,x$Tau),J2=rep(0,x$Tau))
             , B=data.frame(B1=rep(NA,x$Tau),B2=rep(NA,x$Tau))
             , h=data.frame(h1=rep(0,x$Tau),h2=rep(0,x$Tau))
             , shape=rep(0,x$Tau)
             , rate=rep(0,x$Tau)
             , s.diff=rep(0,x$Tau)
             , q.hat=rep(NA,x$Tau)
             , bss=rep(NA,x$Tau)
             , sig=data.frame(sig1=c(rep(0,x$t1.change),rep(1,x$Tau-x$t1.change)),
                              sig2=c(rep(0,x$t2.change),rep(1,x$Tau-x$t2.change)))
             , phi=data.frame(phi1=rep(NA,x$Tau),phi2=rep(NA,x$Tau))
             , Q=data.frame(Q1=rep(0,x$Tau),Q2=rep(0,x$Tau))
             , Pi=data.frame(Pi1=rep(0,x$Tau),Pi2=rep(0,x$Tau))
             , z=list()  
             , s=list() 
             , G=list(G1=list(),G2=list())
             , modelstring=getModelStr(x$params, x$Gdist)
             , t=x$t
             , n.iter = x$n.iter
             , probs = x$probs
             , t1.change = x$t1.change
             , t2.change = x$t2.change
             # , param.inits=x$param.inits
             , Gdist = x$Gdist
  )
  
  t <- l$t
  
  ## CSR strategy costs
  l$phi$phi1 <- sapply(l$sig$sig1, function(s)ifelse(s==1,x$phi1,0))
  l$phi$phi2 <- sapply(l$sig$sig2, function(s)ifelse(s==1,x$phi2,0))
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1,x$gamma2, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  
  for (param in x$params) {
    l$est[[param]] <- data.frame( L99=rep(0,x$Tau),L95=rep(0,x$Tau), mu=rep(0,x$Tau), U95=rep(0,x$Tau),U99=rep(0,x$Tau) )
    l$est[[param]][t, ] <- quantile(rbeta(1e4, x$a1, x$a2), probs = x$probs)
  }
  
  ## Initial values
  l$J$J1[t] <- x$J1.0
  l$J$J2[t] <- x$J2.0
  l$B$B1[t] <- ifelse('B1.0' %in%names(x), x$B1.0, 4 * l$J$J1[t])
  l$B$B2[t] <- ifelse('B2.0' %in%names(x), x$B2.0, 4 * l$J$J2[t])
  l$p$p1[t] <- getPriceFromCsrMarkup(l$sig$sig1[t],l$epsilon,l$gamma$gamma1[t],l$phi$phi1[t], x$c1)   #getPrice(1, x, t)
  l$p$p2[t] <- getPriceFromCsrMarkup(l$sig$sig2[t],l$epsilon,l$gamma$gamma2[t],l$phi$phi2[t], x$c2)  #getPrice(2, x, t)
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, getPsi(l$gamma$gamma2[t],x$Y,l$p$p2[t],l$B$B2[t]), 0)
  l$M[t] <- round(x$db1*l$B$B1[t] + x$db2*l$B$B2[t])
  l$z[[t]] <- rbinom(l$M[t], 1, l$q)  ##  CHANGING TO GROUND TRUTH  q
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
  return(l)
}


##
#
#
##
updateGame <- function(x, l, t)
{
  ## CSR CONTINGENT COSTS
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1, x$gamma1, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  
  ## PERIOD PARAMETERS
  l$p$p1[t] <- getPriceFromCsrMarkup(l$sig$sig1[t],l$epsilon,l$gamma$gamma1[t],l$phi$phi1[t], x$c1)   #getPrice(1, x, t)
  l$p$p2[t] <- getPriceFromCsrMarkup(l$sig$sig2[t],l$epsilon,l$gamma$gamma2[t],l$phi$phi2[t], x$c2)
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
  
  ## CSR-Price-Base-CONTINGENT COSTS
  #cat(sprintf("getPsi: %.10f\n", getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]) ))
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, getPsi(l$gamma$gamma2[t],x$Y,l$p$p2[t],l$B$B2[t]), 0)
  
  # SAMPLE MARKET ATTITUDES
  l$z[[t]] <- rbinom(l$M[t],1, l$q)
  
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
  l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t],dist = l$Gdist)
  l$G$G2[[t]] <- getG(1-l$s[[t]], l$L$L2[t], l$M[t], dist=l$Gdist)
  
  ## LEARN Qhat
  if(x$learn) {
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
learnBayesParams <- function(x,l,t)
{ 
  cat(sprintf('Starting MCMC%s %s iterations\n', ifelse(x$parallel, paste0(' in parallel on ',x$n.cores,' cores'),''),l$n.iter))
  l$s.diff[t] <- mean(l$s[[t]]) - mean(l$s.base)
  s.cumu <- sum(unlist(l$G$G1)) / sum(unlist(c(l$G$G1,l$G$G2)))
  ## Learn customer type by demand groups
  cl <- getDemandClusters(l, t, k='auto', maxcl = 2, all = FALSE)
  l$bss[t] <- getBetweenSS(cl)
  if ( (l$sig$sig1[t] != l$sig$sig2[t]) & (l$bss[t] > x$cl.cutoff) ) {
    plotDemandGroups(l$G$G1[[t]],cl)
    Z <- getZ(x,cl,t)
    h1 <- sum(Z)
    h2 <- length(Z) - h1
    l$q.hat[t] <- h1 / length(Z)
    ## LEARNING params HERE ...
    l$h$h1[t] = ifelse(t==1, x$a1, l$h$h1[t-1]) + h1 * x$ep #downweight
    l$h$h2[t] = ifelse(t==1, x$a2, l$h$h2[t-1]) + h2 * x$ep #downweight
    l$shape[t] = ifelse(t==1,1,l$shape[t-1])  + mean(l$G$G1[t-1])*length(l$G$G1)*x$ep
    l$rate[t] = ifelse(t==1,1,l$rate[t-1])    + length(l$G$G1)*x$ep
  } else {
    ## NOT LEARNING HERE ...
    l$h$h1[t] = ifelse(t==1,x$a1,l$h$h1[t-1])
    l$h$h2[t] = ifelse(t==1,x$a2,l$h$h2[t-1])
    l$shape[t] = ifelse(t==1,1,l$shape[t-1])
    l$rate[t] = ifelse(t==1,1,l$rate[t-1])
  }
  ## -----------------
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
               h2t=l$h$h2[t])
  if ( 'epsilon' %in% x$params) {
    data$shapet <- l$shape[t]
    data$ratet <- l$rate[t] 
  }
  if ( !('epsilon' %in% x$params) ) data$epsilon <- l$epsilon
  if ( !('q' %in% x$params)       ) data$q <- l$q
  
  ##
  # l$param.inits <- list()
  # for (param in x$params) {
  #   l$param.inits[[param]] <- ifelse(param=='q',rbeta(1,l$h$h1[t],l$h$h2[t]),ifelse(param=='epsilon',rgamma(1,1,1),rnorm(1)))
  # }
  
  rngs <- c("base::Super-Duper", "base::Mersenne-Twister",
            "base::Wichmann-Hill", "base::Marsaglia-Multicarry")
  l$param.inits <- list()
  for (i in 1:x$n.cores) {
    l$param.inits[[paste0("inits",i)]] <- list()
    if ('q' %in% x$params) l$param.inits[[paste0("inits",i)]]$q <- rbeta(1,l$h$h1[t],l$h$h2[t])
    if ('epsilon' %in% x$params) l$param.inits[[paste0("inits",i)]]$epsilon <- rgamma(1,1,1)
    l$param.inits[[paste0("inits",i)]]$.RNG.name <- rngs[i]
    l$param.inits[[paste0("inits",i)]]$.RNG.seed <- 12340+i
  }
  
  # inits1 <- list("q"=rbeta(1,l$h$h1[t],l$h$h2[t]), "epsilon"=rgamma(1,1,1), ".RNG.name"=rngs[1], ".RNG.seed"=12341)
  # inits2 <- list("q"=rbeta(1,l$h$h1[t],l$h$h2[t]), "epsilon"=rgamma(1,1,1), ".RNG.name"=rngs[2], ".RNG.seed"=12342)
  # inits3 <- list("q"=rbeta(1,l$h$h1[t],l$h$h2[t]), "epsilon"=rgamma(1,1,1), ".RNG.name"=rngs[3], ".RNG.seed"=12343)
  # inits4 <- list("q"=rbeta(1,l$h$h1[t],l$h$h2[t]), "epsilon"=rgamma(1,1,1), ".RNG.name"=rngs[4], ".RNG.seed"=12344)
  # l$param.inits <- function(){
  #   list(inits1=inits1,inits2=inits2,inits3=inits3,inits4=inits4)
  # }
  
  ### JAGS object if not parallel or MCMC object if parallel
  # if(x$parallel) {
  #   l$sim[[t]] <- runMcmcPar(l,x,t)
  # } else {
  #   l$sim[[t]] <- JAGS$new(l$modelstring, x$params, n.iter.update=x$n.iter)
  #   l$sim[[t]]$run(data, l$param.inits, parallel = x$parallel, period=t) #envir=environment(),
  # }
  #
  l$sim[[t]] <- JAGS$new(l$modelstring, x$params, n.iter=l$n.iter, method=x$method)
  l$sim[[t]]$run(data = data, inits = l$param.inits, parallel = x$parallel,envir = environment(), period=t) #envir=environment(),
  
  
  # if ( !x$parallel & length(x$params) > 1) {
  #   l$sim[[t]]$bivarPlot(x$params[1], x$params[2], chain=1) 
  # }
  if (length(x$params) > 1) {
    bivarPlot(l$sim[[t]]$mcmc, x$params[1], x$params[2], chain=1)
  }
  
  ## estimate and credible intervals
  for (param in x$params) {
    l$est[[param]][t, ] <- getParamEst(l$sim[[t]]$mcmc, param, probs=x$probs, burninProportion = .2)
  }
  
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



#-------------------------------------------------------------------------------------
separateParamCols <- function(df, id.vars, paramCol="param_value") {
  x <- as.character(df[ ,paramCol])
  param.df <- ldply(str_split(str_replace_all(x, "[a-zA-Z]", ""), "_"))
  names(param.df) <- unlist(str_split(str_replace_all(x[1], "[\\d.]", ""), "_"))
  df <- cbind(df[,which(!(names(df) %in% paramCol))], param.df)
  return(df)
}

getBasePlotDf <- function(l.list, id.vars = c('q','epsilon','db','period')) {
  base <- ldply(l.list, function(l) {
    l$B$period <- as.numeric(seq_len(nrow(l$B)))
    l$B$B1share <- l$B$B1 / (l$B$B1 + l$B$B2)
    for(var in names(l$params)) 
      l$B[ ,var] <- l$params[var]
    return(l$B)
  }, .id = "param_value")
  df <- base[,which(!(names(base) %in% "param_value"))]
  base_long <- melt(df,measure.vars = 'B1share',  id.vars = id.vars)
  ##
  df <- subset(base_long, subset = ( !(variable %in% c('B1','B2'))))
  return(df)
}
getSellerPlotDf <- function(l.list, id.vars = c('q','epsilon','db','period')) {
  sell <- ldply(l.list, function(l) {
    l$J$period <- as.numeric(seq_len(nrow(l$J)))
    l$J$J1share <- l$J$J1 / (l$J$J1 + l$J$J2)
    return(l$J)
  }, .id = "param_value")
  df <- separateParamCols(sell, "param_value")
  sell_long <- melt(df, id.vars = id.vars)
  ##
  df <- subset(sell_long, subset = ( !(variable %in% c('J1','J2'))))
  return(df)
}


# ##
# #
# ##
# getModelStr <- function() {
#   "model{
#   for (i in 1:n) {
#   z[i] ~ dbern(q)
#   th1[i] <- p2*(v1 + omega*sig1*z[i])
#   th2[i] <- p1*(v2 + omega*sig2*z[i])
#   s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
#   G[i] ~ dbinom( s[i], L )
#   }
#   q ~ dbeta(h1t,h2t)
# }"
# }

getModelStr <- function(params, dist, vars=c()) {
  for (param in params) {
    vars <- c(vars, switch(param
                           , 'q'='q ~ dbeta(h1t,h2t)'
                           , 'epsilon'='epsilon ~ dgamma(shapet,ratet)')
    )
  }
  Gdist <- switch(dist
                  , 'binomial' = 'dbinom(s[i], L)'
                  , 'poisson' = 'dpois(s[i]*L)' 
  )
  str <- 
    "model{
  for (i in 1:n) {
  z[i] ~ dbern(q)
  th1[i] <- p2*(v1 + omega*sig1*z[i])
  th2[i] <- p1*(v2 + omega*sig2*z[i])
  s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
  G[i] ~ %s
  }
  %s
}"
  return( sprintf(str, Gdist, paste(vars, collapse="\n")) )
  }



##
#
#  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#  'background', 'bgparallel' or 'snow'
###
# x$n.iter = 15000
# 
# system.time(
# rjout <- runjags::run.jags(getModelStr(x$params), monitor=x$params, data=data, 
#                   n.chains=x$n.cores, inits=l$param.inits, 
#                   method="rjparallel",
#                   adapt=ceiling(x$n.iter*(1/15)),
#                   burnin=ceiling(x$n.iter*(4/15)),
#                   sample=ceiling(x$n.iter*(10/15)))
# )
# 
# mcmc <- coda::as.mcmc(mcmc.out)

##
#
##
.sampleMcmcPar <- function(l,x,t, params, n.iter, burninProportion = 0.4)
{ #n.cores = n.chains
  model.file <- "model.bug"
  #jags.inits <- function() {li <- list(); for(i in seq_len(x$n.cores)) li[[i]] <- l$param.inits;  return(li)}
  capture.output(cat(stringr::str_replace_all(l$modelstring,"\n","")), file=model.file)
  samples <- as.mcmc(do.call(jags.parallel,
                             list(data = data, inits=l$param.inits,  
                                  parameters.to.save = x$params,  model.file=model.file,
                                  n.chains=x$n.cores, n.iter=x$n.iter, 
                                  n.burnin=ceiling(x$n.iter*burninProportion))
  )
  )
  return(samples)
}


##
#
##
runMcmcPar <- function(l,x,t, gelmanPlot=F) 
{
  mcmc.output <- .sampleMcmcPar(l, x, t)
  getHeidConvDiag(mcmc.output)
  for (param in x$params) 
    plotMCMCdiagnostics(mcmc.output,t,param)
  if(gelmanPlot)
    gelman.plot(mcmc.output)
  return(list(mcmc = mcmc.output))
}

bivarPlot <- function(mcmc.output, varx, vary, chain=NA, cor.title=TRUE, padx=0, pady=.1) {
  if (length(mcmc.output) > 0) {
    if(is.na(chain) | chain=='all') {
      x <- as.vector(unlist(sapply(mcmc.output,function(x)x[, varx])))
      y <- as.vector(unlist(sapply(mcmc.output,function(x)x[, vary])))
    } else {
      x <- as.vector(mcmc.output[[chain]][, varx])
      y <- as.vector(mcmc.output[[chain]][, vary])          
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
}
