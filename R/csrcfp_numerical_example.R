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
  return(dbinom(xHW, l$N, prob=mu))
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
  if(is.na(rep))  rep <- 1:l$reps
  sigk <- ifelse(k==1, sig1, sig2)
  u1 <- l$u[1]; u2 <- l$u[2];  uk <- l$u[k];
  J1 <- l$J[1]; J2 <- l$J[2];  Jk <- l$J[k]; 
  psi1 <- l$psi.s[1];  psi2 <-l$psi.s[2]; psik <- l$psi.s[k];
  mc1 <- l$c.s[1];  mc2 <- l$c.s[2]; mck <- l$c.s[k]; 
  om<-l$omega; Y<-l$Y; eps<-x$epsilon;
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
  if(is.na(rep))  rep <- 1:l$reps
  sH <- S(l,sig1,sig2,'H',k,rep)
  sU <- S(l,sig1,sig2,'U',k,rep)
  return( sH*(l$z$H[rep] / l$N) + sU*(l$z$U[rep] / l$N) )
}

qstar1.AHAU <- function(l, k=1) {
  Y<-l$Y; N<-l$N; eps<-x$epsilon;   beta <- eps + 1
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
  Y<-l$Y; N<-l$N; eps<-x$epsilon;   beta <- eps + 1
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

init <- function(G, ...) {
  ## update params from inputs
  for (var in names(list(...))) {
    if (var %in% names(G)) {
      G[[var]] <- list(...)[[var]]
    }
  }
  ## Compute derived variables
  z = rbinom(G$reps, size = G$N, G$q)
  G$z <- data.frame(H=z, U=G$N-z)
  G$qstar = qstar(G, k=1:2)
  ## buyers
  G$kstarbar <- data.frame(focal=EKapStar(l, G$qstar$focal, k=1, G$rep),
                           rival=EKapStar(l, G$qstar$rival, k=2, G$rep))
  rownames(G$kstarbar) <- rownames(G$qstar)
  G$sig.b <- sapply(1:ncol(G$kstarbar), function(j){
    sapply(1:nrow(G$kstarbar), function(i){
      ifelse(G$kstarbar[i,j] < G$kap, 1, 0)
    })
  })
  colnames(G$sig.b) <- colnames(G$qstar)
  rownames(G$sig.b) <- rownames(G$qstar)
  ## platform
  G$qhat =  mean(dbeta())
  G$sig.p = data.frame(focal=rep(NA,G$reps), rival=rep(NA,G$reps))
  ## check BPE ?
   # (TODO)
  ## remove init
  G$init = NULL
  return(G)
}


PI <- function(z,sig1,sig2=1, tau=1, N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
               psi.s=.1,  mc.s=1,
               psi.p=.001, mc.p=.01
) {
  s <- S(z,sig1,sig2,k=1)
  beta <- eps + 1
  p1 <- beta * (mc.s + sig1*psi.s)
  marPI <- r - ( (mc.p + psi.p) / p1 )
  out <- sapply(1:tau, function(pd) {
    E <- ifelse(pd==1, E, 0)
    R <- ifelse(pd==1, R, 0)
    if(sig1==sig2) {
      return( pd * (s * (marPI/p1) * N * Y) - E - (R) )
    } else if (z==1) {
      return ( -E )
    } else if (z==0) {
      return( pd * ((marPI/p1) * N * Y) - E )
    }
  })
  return(out)
}

epi <- sapply(qs, function(q){
  csr <- q*PI(1,1,tau=tau) + (1-q)*PI(0,1,tau=tau)
  no <-  q*PI(1,0,tau=tau) + (1-q)*PI(0,0,tau=tau)
  return(c(CSR1=csr, NO1=no))
})
#-----------------------------------------------------------------------------------

## Game definition initializer
G <- list(
  ## G(K,q,N) game parameters
  kap =     .5,   ## buyer's signaling cost
  q =       .2,   ## hedonic proportion
  N =     1000,   ## market size
  reps =   100,   ## game repetitions 
  a = c(.5, .5),  ## prior beliefs about buyer types in market
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
  psi.s = c(0.1, 0.1),         ## seller CSR incremental marginal cost  ##price##
  c.p = c(1, 1),               ## platform marginal cost   
  psi.p = c(0.1, 0.1),         ## platform CSR incremental marginal
  r = c(.1, .1),               ## platform transaction fee rate
  E.s = c(500, 500),           ## seller fixed cost
  R.s = c(100, 100),           ## seller CSR fixed cost
  E.p = c(500, 500),           ## platform fixed cost
  R.p = c(100, 100),           ## platform CSR fixed cost
  ## initialize
  init = function(G, ...){
    G = init(G, ...)
    G$init = NULL
  }
)
 
game <- G$init(q = .2, reps=500)


