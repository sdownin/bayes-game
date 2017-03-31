##-------------------------------------------------------------------------------------
##
##  CSR Signaling Game Numreical Worked Example
##
##-------------------------------------------------------------------------------------
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')
library(lattice)
library(ggplot2)

G <- list(
  ## G(K,q,N) game parameters
  kap = .5, ## buyer's signaling cost
  q =   .2, ## hedonic proportion
  N = 1000, ## market size
  n =  100, ## game repetitions 
  ## strategies
  sig.b = list(W=NA, A=NA),    ## buyer strategies       ## sum to N ## f(z)  
  sig.p = list('1'=0, '2'=1),  ## Platform strategies
  ## decision variables
  z = NA,                      ## r.v. sampled from q,N  
  qstar = NA,                  ## platform's threshold
  qhat = NA,                   ## plaform's estimate of buyers' signals
  ## game variables
  u1=1, u2=1,                  ## baseline (utilitarian) quality of platform
  omega=2,                     ## hedonic quality
  Y=10,                        ## Budget
  J1=80, J2=20,                ## Num Sellers on platform
  epsilon=0.05,                ## Network Effect
  c.s1=1, c.s2=1,              ## seller margina cost by platfrom      ##price##
  psi.s1=0.1, psi.s2=0.1,      ## seller CSR incremental marginal cost ##price##
  c.p=1, c.p2=1,               ## platform marginal cost   
  psi.p1=0.1, psi.p2=0.1,      ## platform CSR incremental marginal
  E.p1 = 500, E.p2=500,        ## fixed cost
  R.p1 = 100, R.p2=100,        ## CSR fixed cost
  ## functions
  init = function(){
    G$z = rbinom(G$n, G$N, G$q)
    return(G)
  }
)

G = G$init()


KapStarBar <- function(zeta,DV,mw,ma) {
  
}

V <- function(sig,z, x) {
  u=1,om=2,Y=10,J=80,eps=0.05,psi=0.1,mc=1
  beta <- eps + 1
  p <- beta * (mc + sig*psi)
  num <- (u + om*sig*z) * Y * J^eps
  return(log2(num/p))
}

qstar <- function(N=1000,Y=100,E=0,R=200, eps=0.05, r=0.1,
                  psi.s=.1,  mc.s=1,
                  psi.p=.001, mc.p=.01) 
{
  beta <- eps + 1
  p1.csr <- beta * (mc.s + psi.s)
  p1.no <- beta * mc.s
  marPI.no <- r - ( (mc.p + psi.p) / p1.no )
  marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
  numer <- N*Y*( (marPI.no/p1.no) - S(0,1,1,1)*(marPI.csr/p1.csr)) + R
  denom <- N*Y*( (marPI.no/p1.no) + (marPI.csr/p1.csr)*(S(1,1,1,1) - S(0,1,1,1)) )
  return(numer / denom)
}

qstarWHAU <- function(N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
                      psi.s=.1,  mc.s=1,
                      psi.p=.001, mc.p=.01
) {
  beta <- eps + 1
  p1.csr <- beta * (mc.s + psi.s)
  p1.no <- beta * mc.s
  marPI.no <- r - ( (mc.p + psi.p) / p1.no )
  marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
  numer <- N*Y*( (marPI.no/p1.no) - (S(1,1,1,1)+S(0,1,1,1))*(marPI.csr/p1.csr)) + 2*R
  denom <- N*Y*( (marPI.no/p1.no) -  S(0,1,1,1)*(marPI.csr/p1.csr) ) + R
  return(numer / denom)
}

## baseline q threshold
png('qstar_exp_profit_by_q.png',height=4.5,width=8, units='in', res=200)
par(mfrow=c(1,2), mar=c(4.2,2,3.5,1))
eps <- c(1.6,  .2);  Rs <- c(20, 2000)
for (i in 1:2) {
  S <- function(z,sig1,sig2,v1=1,v2=1,J1=80,J2=20,omega=2,epsilon=eps[i],k=1,mc=1,psi=0.1
  ) {
    beta <- epsilon + 1
    p1 <- beta * (mc + sig1*psi)
    p2 <- beta * (mc + sig2*psi)
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
  qstar <- function(N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
                    psi.s=.1,  mc.s=1,
                    psi.p=.001, mc.p=.01
  ) {
    beta <- eps + 1
    p1.csr <- beta * (mc.s + psi.s)
    p1.no <- beta * mc.s
    marPI.no <- r - ( (mc.p + psi.p) / p1.no )
    marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
    numer <- N*Y*( (marPI.no/p1.no) - S(0,1,1,1)*(marPI.csr/p1.csr)) + R
    denom <- N*Y*( (marPI.no/p1.no) + (marPI.csr/p1.csr)*(S(1,1,1,1) - S(0,1,1,1)) )
    return(numer / denom)
  }
  qs <- seq(0,1,.01)
  tau <- 1
  
  epi <- sapply(qs, function(q){
    csr <- q*PI(1,1,tau=tau) + (1-q)*PI(0,1,tau=tau)
    no <-  q*PI(1,0,tau=tau) + (1-q)*PI(0,0,tau=tau)
    return(c(CSR1=csr, NO1=no))
  })
  ## plot
  matplot(qs,t(epi), type='l', lty=1:2, lwd=2, col=c('steelblue','darkred'),
          # main=switch(i, '1'=expression('High '*epsilon*' ('*J[1]>J[2]*'), Low'~R),
          #             '2'=expression('Low '*epsilon*', High'~R)),
          main=switch(i, '1'='Low Hedonic Proportion Threshold',
                         '2'='High Hedonic Proportion Threshold'),
          yaxt='n', xlab='Hedonic Proportion (q)',  ylab='')
  title(ylab=expression('Platform Profit'~(pi[1]^p)), line=0, cex.lab=1)
  # abline(h=0, lty=3)
  abline(v=qstar(), lty=4); mtext(text='q*',line = 1, padj = 1 ,at = qstar())
  abline(v=qstarWHAU(), lty=5); mtext(text=expression(q[W^H*A^U]^'*'),line = 1, padj = 1 ,at = qstarWHAU())
  legend(ifelse(i==2,'topright','topright'), title="Strategic Response", 
         legend = c(expression(E[z](pi[1]^p~(CSR)*'|'*q)),
                    expression(E[z](pi[1]^p~(No)*"|"*q))), 
         lty=1:2, col=c('steelblue','darkred'), lwd=2)
}
dev.off()



## posterior q thresholds

png('qstar_exp_profit_by_q.png',height=4.5,width=8, units='in', res=200)
par(mfrow=c(1,2), mar=c(4.2,2,3.5,1))
eps <- c(1.6,  .2);  Rs <- c(20, 2000)
for (i in 1:2) {
  S <- function(z,sig1,sig2,v1=1,v2=1,J1=80,J2=20,omega=2,epsilon=eps[i],k=1,mc=1,psi=0.1
  ) {
    beta <- epsilon + 1
    p1 <- beta * (mc + sig1*psi)
    p2 <- beta * (mc + sig2*psi)
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
  qstar <- function(N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
                    psi.s=.1,  mc.s=1,
                    psi.p=.001, mc.p=.01
  ) {
    beta <- eps + 1
    p1.csr <- beta * (mc.s + psi.s)
    p1.no <- beta * mc.s
    marPI.no <- r - ( (mc.p + psi.p) / p1.no )
    marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
    numer <- N*Y*( (marPI.no/p1.no) - S(0,1,1,1)*(marPI.csr/p1.csr)) + R
    denom <- N*Y*( (marPI.no/p1.no) + (marPI.csr/p1.csr)*(S(1,1,1,1) - S(0,1,1,1)) )
    return(numer / denom)
  }
  qstarWHAU <- function(N=1000,Y=100,E=0,R=Rs[i], eps=0.05, r=0.1,
                        psi.s=.1,  mc.s=1,
                        psi.p=.001, mc.p=.01
  ) {
    beta <- eps + 1
    p1.csr <- beta * (mc.s + psi.s)
    p1.no <- beta * mc.s
    marPI.no <- r - ( (mc.p + psi.p) / p1.no )
    marPI.csr <- r - ( (mc.p + psi.p) / p1.csr )
    numer <- N*Y*( (marPI.no/p1.no) - (S(1,1,1,1)+S(0,1,1,1))*(marPI.csr/p1.csr)) + 2*R
    denom <- N*Y*( (marPI.no/p1.no) -  S(0,1,1,1)*(marPI.csr/p1.csr) ) + R
    return(numer / denom)
  }
  qs <- seq(0,1,.01)
  tau <- 1
  
  epi <- sapply(qs, function(q){
    csr <- q*PI(1,1,tau=tau) + (1-q)*PI(0,1,tau=tau)
    no <-  q*PI(1,0,tau=tau) + (1-q)*PI(0,0,tau=tau)
    return(c(CSR1=csr, NO1=no))
  })
  ## plot
  matplot(qs,t(epi), type='l', lty=1:2, lwd=2, col=c('steelblue','darkred'),
          # main=switch(i, '1'=expression('High '*epsilon*' ('*J[1]>J[2]*'), Low'~R),
          #             '2'=expression('Low '*epsilon*', High'~R)),
          main=switch(i, '1'='Low Hedonic Proportion Threshold',
                      '2'='High Hedonic Proportion Threshold'),
          yaxt='n', xlab='Hedonic Proportion (q)',  ylab='')
  title(ylab=expression('Platform Profit'~(pi[1]^p)), line=0, cex.lab=1)
  # abline(h=0, lty=3)
  abline(v=qstar(), lty=4); mtext(text='q*',line = 1, padj = 1 ,at = qstar())
  abline(v=qstarWHAU(), lty=5); mtext(text=expression(q[W^H*A^U]^'*'),line = 1, padj = 1 ,at = qstarWHAU())
  legend(ifelse(i==2,'topright','topright'), title="Strategic Response", 
         legend = c(expression(E[z](pi[1]^p~(CSR)*'|'*q)),
                    expression(E[z](pi[1]^p~(No)*"|"*q))), 
         lty=1:2, col=c('steelblue','darkred'), lwd=2)
}
dev.off()


q <- .6
tau <- 10
df <- data.frame(
  CSR=q*PI(1,1,tau=tau) + (1-q)*PI(0,1,tau=tau),
  NO=q*PI(1,0,tau=tau) + (1-q)*PI(0,0,tau=tau)
)
matplot(df, type='l', lwd=2)
abline(h=0, lty=4)

PI(1,1)
PI(0,1)
PI(1,0)
PI(0,0)

## BUYER VALUE
# V(1,1) - kap
V(1,1)
# V(1,0) - kap
V(1,0) 
# V(0,1) - kap
V(0,1)
# V(0,0) - kap
V(0,0) 

df <- data.frame(
  share=c(S(1,1,0,k=1), S(1,0,0,k=1), S(0,1,0,k=1), S(0,0,0,k=1),
          S(1,1,1,k=1), S(1,0,1,k=1), S(0,1,1,k=1), S(0,0,1,k=1),
          S(1,0,1,k=2), S(1,0,0,k=2), S(0,0,1,k=2), S(0,0,0,k=2),
          S(1,1,1,k=2), S(1,1,0,k=2), S(0,1,1,k=2), S(0,1,0,k=2))
)
df$z <- factor(rep(c('H','H','U','U'),4))
df$strategy <- factor(c(rep(c('1 = CSR','1 = NO'),4),rep('1 = NO',4),rep('1 = CSR',4)))
df$rival <- factor(c(rep('2 = NO',4),rep('2 = CSR',4),rep(c('2 = CSR','2 = NO'),4)))
df$platform <- factor(c(rep('1',8),rep('2',8)))
df$z.strat <- sapply(1:nrow(df), function(x)paste(df$z[x], df$strategy[x], sep="_"))

sapply(levels(df$strategy), function(i) {
  sapply(levels(df$rival), function(j) {
    sub <- df[df$strategy==i & df$rival==j, ]
    return(c(Es1=mean(sub$share[sub$platform=='1']), 
             Es2=mean(sub$share[sub$platform=='2'])))
  })
})

# lattice::dotplot(share ~ z.strat | platform, groups=rival, data=df, 
#                  auto.key=T, pch=16, ylim=c(0,1))

ggplot(aes(x=z, y=share, fill=platform), data=df) + 
  geom_bar(aes(fill=platform), position='dodge', stat='identity') + 
  facet_grid(strategy ~ rival) + 
  ylim(0,1) + theme_bw()




PSI <- function(kap=1,om=2,u=1,c=1) {
  ratio <- (u + om) / (u*exp(kap))
  return( c * (ratio - 1) ) 
}

kaps <- seq(0,2,.01)
psis <- PSI(kap=kaps)

par(mar=c(4.2,4.2,3.5,1))
plot(kaps, psis, type='l', lwd=2, 
     xlab='Buyer"s Signaling Cost', 
     ylab='Seller"s Marginal CSR Cost')
polygon(x=c(0,kaps),y=c(min(psis),psis), angle=90,col=rgb(.2,.2,.2,.3))
abline(h=0)

KAP <- function(om,psi,u=2,c=1) {
  value <-(u+om) / u
  cost <- c / (c + psi)
  return(log(value*cost))
}

oms <- seq(0,2,.005)
psis <- seq(0,2,.005)

VecKAP <- Vectorize(FUN = KAP)

xx <- outer(oms, psis, VecKAP)

image(xx); contour(xx, add=T)



PrW <- function(Y=100, u=1, mc.s=1, J1=80, eps=1.05) {
  p.no <- (1+eps) * mc.s
  ratio <- p.no / (u * Y * J1^eps)
  return( -log(ratio) )
}


PrW <- function(kap=1, Y=100, u=1, mc.s=1, psi.s=0.1, J1=80, eps=1.05, om=2) {
  p.no <- (1+eps) * mc.s
  p.csr <- (1+eps) * (mc.s + psi.s)
  value <- (u+om) / u
  numer <- log(value * (p.no/p.csr) ) - kap
  denom <- log( ((u+om)*Y*J1^eps)/p.csr ) - kap
  return( numer/denom )
}

MU <- function(kap=0.5, Y=100, u=1, mc.s=1, psi.s=0.1, J1=80, eps=1.05, om=2) {
  p.no <- (1+eps) * mc.s
  p.csr <- (1+eps) * (mc.s + psi.s)
  value <- (u+om) / u
  zeta <- (u+om) / (p.csr + exp(kap))
  numer <- log((p.no/u) * zeta)
  denom <- log(Y*J1^eps * zeta )
  return( numer/denom )
}


cdf <- sapply(0:10, function(k){
  x <- 1:k
  y <- sum(choose(N,x) * beta(x+a1, N-x+a2)/beta(a1,a2) )
  return(y)
})
plot(0:10,cdf)

dbinom(1, 1, prob = MU() )

N <- 1000
a1 <- a2 <- .5
qstarN <- qstar() * (N+a1+a2) - a1
1 - sum(dbinom(1:qstarN,N,MU()))
 
MU <- function(kap=.5, u=1, om=2, c.s=1, psi.s=0.1) {
  value <- (u+om)/u
  cost <- c.s / (c.s+psi.s)
  return( 1- (kap / log(value * cost)) )
}

kaps <- seq(0,3,.01)
plot(MU(kap=kaps,om=5) ~ kaps, ylim=c(0,1), 
     ylab=expression(P(sigma[i]^b==W*'|'*kappa)),
     xlab=expression('Signal Cost'~(kappa)),
     type='l')





kapstar <- function(dV, )