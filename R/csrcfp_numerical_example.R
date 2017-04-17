##-------------------------------------------------------------------------------------
##
##  CSR Signaling Game Numreical Worked Example
##
##-------------------------------------------------------------------------------------
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')
library(lattice)
library(ggplot2)

#-----------------------------------------------------------------------------------
xHW <- function(qstar, l) { 
  val <- qstar*(l$N + l$a[1] + l$a[2]) - l$a[1] + 1
  return(round(val))
}

MU <- function(l, k=1) {
  c.s <- l$c.s[k];   psi.s <- l$psi.s[k]; u <- l$u[k]; J <- l$J[k];
  om <- l$omega; Y <- l$Y; kap <- l$kap;  eps <- l$epsilon; beta <- eps+1;
  ekap <- 1/exp(kap)
  top <- ekap * (c.s/(c.s+psi.s)) * ((u+om)/u)
  bottom <- ekap * (Y*J^eps * (u+om)) / (beta*(c.s+psi.s))
  return(log(top)/log(bottom))
}

V <- function(l, k=1, z='H', sig=0, rep=1) {
  u <- l$u[k]; J <- l$J[k]; psi <- l$psi.s[k]; mc <- l$c.s[k];
  om<-l$omega; Y<-l$Y; eps<-l$epsilon;  
  zval <- ifelse(z=='H',1,0)
  beta <- eps + 1
  p <- beta * (mc + sig*psi)
  num <- (u + om*sig*zval) * Y * J^eps
  return(log(num/p))
}

mdiff <- function(xHW, N, mu, k=1) {
  return(dbinom(xHW, N, prob=mu))
}

Vdiff <- function(l, k=1, rep=1) {
  return(V(l,k,z='H',sig=1,rep) - V(l,k,z='U',sig=1,rep))
}

EKapStar <- function(l, qstar, k=1, rep=1) {
  .xHW <- xHW(qstar, l)
  .mu <- MU(l, k)
  .mdiff <- mdiff(.xHW, l$N, .mu)
  .Vdiff <- Vdiff(l, k, rep)
  return(l$q * .Vdiff * .mdiff)
}

S <- function(l, sig1, sig2, z='H', k=1, rep=1) {
  if(any(is.na(rep)))  rep <- 1:l$reps
  sigk <- ifelse(k==1, sig1, sig2)
  u1 <- l$u[1]; u2 <- l$u[2];  uk <- l$u[k];
  J1 <- l$J[1]; J2 <- l$J[2];  Jk <- l$J[k]; 
  psi1 <- l$psi.s[1];  psi2 <-l$psi.s[2]; psik <- l$psi.s[k];
  mc1 <- l$c.s[1];  mc2 <- l$c.s[2]; mck <- l$c.s[k]; 
  om<-l$omega; Y<-l$Y; eps<-l$epsilon;
  beta <- eps + 1
  ## prices
  p1 <- beta * (mc1 + sig1*psi1)
  p2 <- beta * (mc2 + sig2*psi2)
  pk <- ifelse(k==1, p2, p1)
  ##
  zval <- ifelse(z=='H',1,0)
  ## thetas
  th1 <- p2 * (u1 + om*sig1*zval)
  th2 <- p1 * (u2 + om*sig2*zval)
  thk <- pk * (uk + om*sigk*zval)
  return( (thk * Jk^eps) / ((th1 * J1^eps)+(th2 * J2^eps)) )
}

ES <- function(l, sig1, sig2, k=1, rep=1) {
  if(any(is.na(rep)))  rep <- 1:l$reps
  sH <- S(l,sig1,sig2,'H',k,rep)
  sU <- S(l,sig1,sig2,'U',k,rep)
  return( sH*(l$z$H[rep] / l$N) + sU*(l$z$U[rep] / l$N) )
}

qstar1.AHAU <- function(l, k=1) {
  Y<-l$Y; N<-l$N; eps<-l$epsilon;   beta <- eps + 1
  r <- l$r[k]
  psi.s <- l$psi.s[k]; mc.s <- l$c.s[k]
  psi.p <- l$psi.p[k]; mc.p <- l$c.p[k]
  E <- l$E.p[k]; R <- l$R.p[k]
  ##
  p1.csr <- beta * (mc.s + psi.s)
  p1.no <- beta * mc.s
  marPI.no <- r - ( (mc.p + psi.p) / p1.no )
  marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
  numer <- N*Y*( (marPI.no/p1.no) - S(l,1,1,'U',k,rep)*(marPI.csr/p1.csr)) + R
  denom <- N*Y*( (marPI.no/p1.no) + (marPI.csr/p1.csr)*(S(l,1,1,'H',k,rep) - S(l,1,1,'U',k,rep)) )
  return(numer / denom)
}

qstar2.WHAU <- function(l, k=1) {
  Y<-l$Y; N<-l$N; eps<-l$epsilon;   beta <- eps + 1
  r <- l$r[k]
  psi.s <- l$psi.s[k]; mc.s <- l$c.s[k]
  psi.p <- l$psi.p[k]; mc.p <- l$c.p[k]
  E <- l$E.p[k]; R <- l$R.p[k]
  ##
  beta <- eps + 1
  p1.csr <- beta * (mc.s + psi.s)
  p1.no <- beta * mc.s
  marPI.no <- r - ( (mc.p + psi.p) / p1.no )
  marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
  numer <- N*Y*( (marPI.no/p1.no) - (S(l,1,1,'H',k,rep)+S(l,1,1,'U',k,rep))*(marPI.csr/p1.csr)) + 2*R
  denom <- N*Y*( (marPI.no/p1.no) -  S(l,1,1,'U',k,rep)*(marPI.csr/p1.csr) ) + R
  return(numer / denom)
}

qstar3.AHWU <- function(l, k=1) {
  return(sapply(k, function(k_i)0))
}
qstar4.WHWU <- function(l, k=1) {
  return(sapply(k, function(k_i)-Inf))
}

qstar <- function(l, k=1) {
  df <- data.frame(sapply(k, function(k_i) {
    .q1 <- qstar1.AHAU(l, k_i)
    .q2 <- qstar2.WHAU(l, k_i)
    .q3 <- qstar3.AHWU(l, k_i)
    .q4 <- qstar4.WHWU(l, k_i)
    return(c(AHAU=.q1, WHAU=.q2, AHWU=.q3, WHWU=.q4))
  }))
  names(df) <- c('focal','rival')[k]
  return(df)
}

#-----------------------------------------
# Strategy Decisions
#-----------------------------------------
init.buyers <- function(G, ...) {
  ## threshold: buyer signaling cost-to-benefit
  G$kstarbar <- data.frame(focal=EKapStar(G, G$qstar$focal, k=1, G$rep),
                           rival=EKapStar(G, G$qstar$rival, k=2, G$rep))
  rownames(G$kstarbar) <- rownames(G$qstar)
  G$sig.b <- sapply(1:ncol(G$kstarbar), function(j){
      sapply(1:nrow(G$kstarbar), function(i){
        return(ifelse(G$kap < G$kstarbar[i,j], 1, 0))
      })
    })
  colnames(G$sig.b) <- colnames(G$qstar); rownames(G$sig.b) <- rownames(G$qstar)
  ## number of buyers voting if ALL H buyers actually vote
  G$sig.b.sum <- lapply(G$z$H, function(x) G$sig.b * x)
  .mu <- MU(G, k=1:2)
  ## number of buyers voting if H buyers vote with prob. = mu_k
  G$sig.b.sum.r <- lapply(G$z$H, function(x) {
    tmp <- G$sig.b
    tmp[,1] <- G$sig.b[,1] * rbinom(1, x, MU(G, k=1))
    tmp[,2] <- G$sig.b[,2] * rbinom(1, x, MU(G, k=2))
    return(tmp)
  })
  return(G)
}

init.platform <- function(G, ...) {
  ## threshold: firm CSR cost-to-benefit
  ###
  # sig.b.sum or sig.b.sum.r ?
  ##
  buyerVoteList <- G$sig.b.sum.r  ## G$sig.b.sum   ### ??
  G$qhat <-  lapply(buyerVoteList, function(sig.b.sum.r){
    tmp <- (sig.b.sum.r + G$a[1]) / (G$N + G$a[1] + G$a[2])
    colnames(tmp) <- colnames(G$qstar); rownames(tmp) <- rownames(G$qstar)
    return(tmp)
  })
  ##
  G$qhat.r <-  lapply(buyerVoteList, function(sig.b.sum.r){
    tmp <- (sig.b.sum.r + G$a[1]) / (G$N + G$a[1] + G$a[2])
    colnames(tmp) <- colnames(G$qstar); rownames(tmp) <- rownames(G$qstar)
    return(tmp)
  })
  ##
  G$qhat.mean <-  as.data.frame(sapply(1:ncol(G$kstarbar), function(j){
    sapply(1:nrow(G$kstarbar), function(i){
      return(mean(sapply(seq_along(G$qhat), function(k){
        G$qhat[[k]][i,j]
      }), na.rm = T))
    })
  }))
  colnames(G$qhat.mean) <- colnames(G$qstar); rownames(G$qhat.mean) <- rownames(G$qstar)
  ##
  G$qhat.r.mean <-  as.data.frame(sapply(1:ncol(G$kstarbar), function(j){
    sapply(1:nrow(G$kstarbar), function(i){
      return(mean(sapply(seq_along(G$qhat.r), function(k){
        G$qhat[[k]][i,j]
      }), na.rm = T))
    })
  }))
  colnames(G$qhat.r.mean) <- colnames(G$qstar); rownames(G$qhat.r.mean) <- rownames(G$qstar)
  ## firm CSR strategy
  G$sig.p <- lapply(seq_along(G$qhat), function(i){
    qhat <- G$qhat[[i]]
    return( (qhat > G$qstar)*1 )  ## *1 changes bool to numeric while keeping data.frame structure
  })
  ## mean of firm CSR strategy over repeated simulations
  G$sig.p.mean <- as.data.frame(sapply(1:ncol(G$kstarbar), function(j){
    sapply(1:nrow(G$kstarbar), function(i){
      return(mean(sapply(seq_along(G$sig.p), function(k){
        G$sig.p[[k]][i,j]
      }), na.rm = T))
    })
  }))
  colnames(G$sig.p.mean) <- colnames(G$qstar); rownames(G$sig.p.mean) <- rownames(G$qstar)

  return(G)
}


# PI <- function(z,sig1,sig2=1, tau=1, N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
#                psi.s=.1,  mc.s=1,
#                psi.p=.001, mc.p=.01
# ) {
#   s <- S(z,sig1,sig2,k=1)
#   beta <- eps + 1
#   p1 <- beta * (mc.s + sig1*psi.s)
#   marPI <- r - ( (mc.p + psi.p) / p1 )
#   out <- sapply(1:tau, function(pd) {
#     E <- ifelse(pd==1, E, 0)
#     R <- ifelse(pd==1, R, 0)
#     if(sig1==sig2) {
#       return( pd * (s * (marPI/p1) * N * Y) - E - (R) )
#     } else if (z==1) {
#       return ( -E )
#     } else if (z==0) {
#       return( pd * ((marPI/p1) * N * Y) - E )
#     }
#   })
#   return(out)
# }
# 
# epi <- sapply(qs, function(q){
#   csr <- q*PI(1,1,tau=tau) + (1-q)*PI(0,1,tau=tau)
#   no <-  q*PI(1,0,tau=tau) + (1-q)*PI(0,0,tau=tau)
#   return(c(CSR1=csr, NO1=no))
# })
#-----------------------------------------------------------------------------------




#-------------------------------------------------------------
# Single-Period Game
#-------------------------------------------------------------

## Game definition initializer
G <- list(
  ## G(K,q,N) game parameters
  kap =     .001,   ## buyer's signaling cost
  q =       .2,   ## hedonic proportion
  N =     10,   ## market size
  reps =   100,   ## game repetitions 
  a = c(.5, .5),  ## prior beliefs about buyer types in market
  seed = 1111,
  ## strategies
  sig.b = c(),    ## buyer strategies       ## sum to N ## f(z)  
  sig.p = c(),    ## Platform strategies
  ## decision variables
  z = c(),                     ## r.v. sampled from q,N  
  kstarbar = c(),              ## expected buyer's signaling cost threshold
  qstar =c(),                  ## platform's threshold
  qhat = c(),                  ## plaform's estimate of buyers' signals
  ## game variables
  omega=2,                     ## hedonic quality
  Y=10,                        ## Budget
  epsilon=1.05,                ## Network Effect
  u = c(1, 1),                 ## baseline (utilitarian) quality of platform
  J = c(80, 20),               ## Num Sellers on platform
  c.s = c(1, 1),               ## seller marginal cost                  ##price##
  psi.s = c(0.01, 0.01),         ## seller CSR incremental marginal cost  ##price##
  c.p = c(1, 1),               ## platform marginal cost   
  psi.p = c(0.01, 0.01),         ## platform CSR incremental marginal
  r = c(.1, .1),               ## platform transaction fee rate
  E.s = c(500, 500),           ## seller fixed cost
  R.s = c(100, 100),           ## seller CSR fixed cost
  E.p = c(500, 500),           ## platform fixed cost
  R.p = c(100, 100),           ## platform CSR fixed cost
  ## initialize
  init = function(...) {
    for (var in names(list(...))) {  ## update params from inputs
      if (var %in% names(G)) {
        G[[var]] <- list(...)[[var]]
      }
    }
    set.seed(G$seed)
    G$qstar = qstar(G, k=1:2)  ## doesn't change per repetition
    z = rbinom(G$reps, size = G$N, G$q)
    G$z <- data.frame(H=z, U=G$N-z)
    G <- init.buyers(G)
    G <- init.platform(G)
    ## check BPE ? 
      ## TODO
    G$init = NULL
    return(G)
  }
)
 





library(lattice)
library(grid)
library(gridExtra)
library(cowplot)
##------------------------------------------------
##  Simulation
##------------------------------------------------
.reps <- 500
.q <- .25
#---------------------------------------
pl <- list()
Ns <- c(10,100, 1000, 10000)
for (p_i in seq_along(Ns)) {
  N <- Ns[p_i]
  g4 <- G$init(q = .q, N=N, reps=.reps)
  ## market distribution
  cH <- plyr::count(g4$z$H); cH$Z <- 'H'; cH$N <- factor(N)
  cU <- plyr::count(g4$z$U); cU$Z <- 'U'; cU$N <- factor(N)
  dfc4 <- rbind(cH,cU)
  dfc4$freq <- round( dfc4$freq / g4$reps, 5)
  ## qstar
  cat(sprintf('\nN = %s\n      --------qstar-----------firm--------buyer----\n',N));
  AU <- 1:4
  print(cbind(g4$qstar[AU,], g4$sig.p.mean[AU,], g4$sig.b[AU,]))
  ## BARPLOT
  pl[[p_i]] <- ggplot(aes(x, freq, fill=Z), data=dfc4)  +
    geom_bar(stat='identity') +
     # geom_text(aes(label=freq), vjust=1.6, color="white",
     #          position = position_dodge(0.9), size=3/log10(g4$N)) +
    geom_vline(xintercept=g4$qstar$focal*g4$N, lty=1:4) +
    geom_vline(xintercept=g4$qhat.r.mean$focal*g4$N, lty=1, colour='green') +
    geom_vline(xintercept=g4$qhat.mean$focal*g4$N, lty=1, colour='green') +
    ylab(sprintf('Mean Sample Mass', g4$reps)) + 
    xlim(-1,g4$N+1) +
    xlab('Sample Count of Buyers by Type') +
    ggtitle(sprintf('N = %s',g4$N)) +
    theme_bw() 
}

########################################################
plot_grid(pl[[1]],pl[[2]], pl[[3]], pl[[4]], 
          labels=c('A','B', 'C','D'), ncol=2, nrow=2)
#######################################################
# 
# .reps <- 500
# .q <- .18
# #---------------------------------------
# pl <- list()
# Ns <- c(10,100, 1000, 10000)
# for (p_i in seq_along(Ns)) {
#   N <- Ns[p_i]
#   g4 <- G$init(q = .q, N=N, reps=.reps)
#   ## market distribution
#   cH <- plyr::count(g4$z$H); cH$Z <- 'H'; cH$N <- factor(N)
#   cU <- plyr::count(g4$z$U); cU$Z <- 'U'; cU$N <- factor(N)
#   dfc4 <- rbind(cH,cU)
#   dfc4$freq <- round( dfc4$freq / g4$reps, 2)
#   ## qstar
#   cat(sprintf('\nN = %s\n      --------qstar-----------firm--------buyer----\n',N));
#   print(cbind(g4$qstar,g4$sig.p[[1]], g4$sig.b))
#   ## BARPLOT
#   pl[[p_i]] <- ggplot(aes(x, fill=Z), data=dfc4)  +
#     geom_histogram(bins=50) +
#     geom_vline(xintercept=g4$qstar$focal*g4$N, lty=1:4, lwd=1.2) +
#     ylab(sprintf('Mean Sample Mass', g4$reps)) + 
#     xlim(-1,g4$N+1) +
#     xlab('Buyer Distribution') +
#     ggtitle(sprintf('N = %s',g4$N)) +
#     theme_bw() 
# }
# 
# ########################################################
# plot_grid(pl[[1]],pl[[2]], pl[[3]], pl[[4]], 
#           labels=c('A','B', 'C','D'), ncol=2, nrow=2)
# #######################################################


# barchart(freq ~ factor(x), groups= Z, data=dfc, auto.key=T)
# 
# cHU <- merge(cH,cU,by='x',all=T)
# barplot(as.matrix(cHU[,-1]), names.arg = cHU$x,  beside = T)
# hist(g1$z$H, xlim=c(0,g1$N), col=rgb(.2,.2,.8,.4)); hist(g1$z$U, add=T, col=rgb(.8,.2,.2,.4))

## decision thresholds
g1$qstar
g1$kstarbar

## strategies
g1$sig.b
g1$sig.p[[1]]

















#----------------------------------------------------------
#  Simulation Individual Charts
#-----------------------------------------------------------

# g1 <- G$init(q = .q, N=10, reps=.reps)
# ## market distribution
# cH <- plyr::count(g1$z$H); cH$Z <- 'H'; cH$N <- factor(N)
# cU <- plyr::count(g1$z$U); cU$Z <- 'U'; cU$N <- factor(N)
# dfc1 <- rbind(cH,cU)
# dfc1$freq <- round( dfc1$freq / g1$reps, 2)
# ## BARPLOT
# p1 <- ggplot(aes(x, freq, fill=Z), data=dfc1)  +
#   geom_bar(stat='identity', position=position_dodge()) +
#    geom_text(aes(label=freq), vjust=1.6, color="white",
#             position = position_dodge(0.9), size=3/log10(g1$N)) + 
#   ylab(sprintf('Mean Sample Mass', g1$reps)) +
#   xlab('Buyer Distribution') +
#   ggtitle(sprintf('N = %s',g1$N)) + 
#   theme_bw()
# #----------------------------------
# g2 <- G$init(q = .q, N=20, reps=.reps)
# ## market distribution
# cH <- plyr::count(g2$z$H); cH$Z <- 'H'; cH$N <- factor(N)
# cU <- plyr::count(g2$z$U); cU$Z <- 'U'; cU$N <- factor(N)
# dfc2 <- rbind(cH,cU)
# dfc2$freq <- round( dfc2$freq / g2$reps, 2)
# ## BARPLOT
# p2 <- ggplot(aes(x, freq, fill=Z), data=dfc2)  +
#   geom_bar(stat='identity', position=position_dodge()) +
#    geom_text(aes(label=freq), vjust=1.6, color="white",
#             position = position_dodge(0.9), size=3/log10(g2$N)) + 
#   ylab(sprintf('Mean Sample Mass', g2$reps)) +
#   xlab('Buyer Distribution') +
#   ggtitle(sprintf('N = %s',g2$N)) + 
#   theme_bw()
# #----------------------------------
# g3 <- G$init(q = .q, N=40, reps=.reps)
# ## market distribution
# cH <- plyr::count(g3$z$H); cH$Z <- 'H'; cH$N <- factor(N)
# cU <- plyr::count(g3$z$U); cU$Z <- 'U'; cU$N <- factor(N)
# dfc3 <- rbind(cH,cU)
# dfc3$freq <- round( dfc3$freq / g3$reps, 2)
# ## BARPLOT
# p3 <- ggplot(aes(x, freq, fill=Z), data=dfc3)  +
#   geom_bar(stat='identity', position=position_dodge()) +
#    geom_text(aes(label=freq), vjust=1.6, color="white",
#             position = position_dodge(0.9), size=3/log10(g3$N)) + 
#   ylab(sprintf('Mean Sample Mass', g3$reps)) +
#   xlab('Buyer Distribution') +
#   ggtitle(sprintf('N = %s',g3$N)) + 
#   theme_bw()
# #----------------------------------
# g4 <- G$init(q = .q, N=80, reps=.reps)
# ## market distribution
# cH <- plyr::count(g4$z$H); cH$Z <- 'H'; cH$N <- factor(N)
# cU <- plyr::count(g4$z$U); cU$Z <- 'U'; cU$N <- factor(N)
# dfc4 <- rbind(cH,cU)
# dfc4$freq <- round( dfc4$freq / g4$reps, 2)
# ## BARPLOT
# p4 <- ggplot(aes(x, freq, fill=Z), data=dfc4)  +
#   geom_bar(stat='identity', position=position_dodge()) +
#    geom_text(aes(label=freq), vjust=1.6, color="white",
#             position = position_dodge(0.9), size=3/log10(g4$N)) + 
#   ylab(sprintf('Mean Sample Mass', g4$reps)) +
#   xlab('Buyer Distribution') +
#   ggtitle(sprintf('N = %s',g4$N)) + 
#   theme_bw()