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
      th1 <- x$p2 * (x$v1 + x$w * x$sig1 * z)
      th2 <- x$p1 * (x$v2 + x$w * x$sig2 * z)
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

share <- function(p1,p2,v1,v2,sig1,sig2,J1,J2,w,z,rho,k=1) 
{
  out <- sapply(z, function(z_i){
    th1 <- p2 * (v1 + w * sig1 * z_i)
    th2 <- p1 * (v2 + w * sig2 * z_i)
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
                  sig1 = 1, sig2 = 0, w = 1.5 , rho = 1, z = NA)
{
  return(list(v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,sig1=sig1,sig2=sig2,w=w,rho=rho,z=z))
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

#------------------ SIMULATE FROM GROUND TRUTH q=.2-----------------------------


df.sim <- getSimulation(q=.5,Y=10000, N=1000)

# ## OVERLAPPING HISTOGRAMS
# z1 <- ifelse(sum(df.sim$G1[df.sim$z==1])>sum(df.sim$G1[df.sim$z==0]), 1, 0)
# z2 <- 1-z1
# hist(df.sim$G1[df.sim$z==z1],col=rgb(1,.1,.1,.5))
# hist(df.sim$G1[df.sim$z==z2],col=rgb(.1,.1,1,.5),add=TRUE)

# ## GGPLOT
# df.sim$z <- as.factor(df.sim$z)
# p1 <- ggplot(aes(x=G1), data=df.sim) +  geom_density(aes(group=z, fill=z),alpha=.5, bw=.5)
# p1
# p2 <- ggplot(aes(x=G1), data=df.sim) +  geom_histogram(aes(group=z, fill=z),alpha=.5, bins=10)
# p2


## COMPARE DIFFERENT q effect
df <- data.frame()
for (q in c(0.05, .25,.75,.95)) {
  set.seed(111)
  df.sim <- getSimulation(q=q, N=5000, Y=20000)
  df.sim$q <- as.factor(paste0('q = ',q))
  df <- rbind(df,df.sim)
}
df$z <- as.factor(df$z)

## FACET WRAP
## GGPLOT

# p1 <- ggplot(aes(x=G1), data=df) +
#   geom_density(aes(group=z, fill=z),alpha=.5, bw=.8) +
#   facet_wrap( ~ q, nrow=2 )
# p1
#--------------------------------------------
## GGPLOT CURRENTLY BROKEN ???
# p2 <- ggplot(aes(x=G1, data=df)) +  
#   geom_histogram(aes(group=z, colour=z, fill=z),alpha=1, bins=15) +
#   facet_wrap( ~ q, nrow=2 ) + 
#   xlab('Quantity Demanded from Platform 1') + 
#   scale_fill_manual(values=c("#003399", "#99bbff")) + 
#   scale_colour_manual(values=c("#001133", "#0044cc")) +
#   theme_bw()
# p2
# ggsave('demand_by_q_4facet.png',width = 6, height = 4.5, units = 'in', dpi=250)
#-------------------------------


## COMPARE DIFFERENT q effect
## Run 100 simulations each at 4 value s of q
runs <- 60
qvec <- seq(0,1,.1)
df.comb <- data.frame()
set.seed(111)
for (run in 1:runs) {
    df <- data.frame() 
    for (q in qvec) {
        df.sim <- getSimulation(q=q, N=60, Y=200, setSeed = FALSE)
        df.sim$q <- as.factor(q)
        df <- rbind(df,df.sim)
    }
    df$z <- as.factor(paste0('z=',df$z))
    df.2 <- rbind(
        data.frame(Demand=as.integer(df$G1),
                   z=df$z,
                   q=df$q,
                   Platform='Platform=1'),
        data.frame(Demand=as.integer(df$G2),
                   z=df$z,
                   q=df$q,
                   Platform='Platform=2')
    )
    df.counts <- ddply(.data = df.2, .variables = .(z,q,Platform),summarise,counts=sum(Demand))
    df.comb <- rbind(df.comb,df.counts)
}

## BW PLOT
png('csr_demand_q_z_bwplot.png', height=6,width=8, units='in', res=300)
bwplot(counts ~ q | z + Platform , 
       data=df.comb,
       xlab='Proportion of Hedonic Demand (q)',
       col='black',
       par.settings = list(strip.background=list(col="lightgrey"),
                           box.umbrella=list(col='black'),
                           box.rectangle=list(col='black'),
                           box.dot=list(col='black')
                           )
       )
dev.off()

##-------------- PLOT SIMULATION ---------------------------------------
## Run 100 simulations each at 4 value s of q
runs <- 100
N <- 500
Y <- 1000
qvec <- seq(0,1,.025)
df.comb <- data.frame()
li.comb <- list()
set.seed(111)
for (run in 1:runs) {
  df <- data.frame() 
  for (q in qvec) {
    df.sim <- getSimulation(q=q, N=N, Y=Y, setSeed = FALSE)
    df.sim$q <- q
    df <- rbind(df,df.sim)
  }
  df$z <- as.factor(paste0(ifelse(df$z==1,'Hedonic ','Utilitarian '),'(z=',df$z,')'))
  df.2 <- rbind(
    data.frame(Demand=as.integer(df$G1),
               z=df$z,
               q=df$q,
               Platform='1 (CSR)'),
    data.frame(Demand=as.integer(df$G2),
               z=df$z,
               q=df$q,
               Platform='2 (NO CSR)')
  )
  df.counts <- ddply(.data = df.2, .variables = .(z,q,Platform),summarise,counts=sum(Demand))
  li.comb[[run]] <- df.counts
  df.comb <- rbind(df.comb,df.counts)
}
# values=c(rgb(.8,.8,.5,.7), rgb(.1,.1,.5,.7))

# GET RIBBON 95% CREDIBLE INTERVALS
df.ci <- data.frame(z=NA, q=NA, Platform=NA, L99=NA, mu=NA, U99=NA)
counter <- 1
for (i in levels(df.counts$z)) {
  for (j in unique(df.counts$q)) {
    for (k in levels(df.counts$Platform)) {
      x <- df.comb[which(df.comb$z==i & df.comb$q==j & df.comb$Platform==k), 'counts']
      df.ci[counter,c('z','q','Platform')] <- c(i,j,k)
      df.ci[counter,c('L99','mu','U99')] <- quantile(na.omit(x), probs = c(.005, .5, .995))
      counter <- counter + 1
    }
  }
}

p <- ggplot(aes(x=as.numeric(q),y=mu), data=na.omit(df.ci)) +
  # geom_smooth(aes(x=q, y=counts, colour=z), data=df.comb) +
  geom_point(aes(colour=Platform, shape=Platform)) +
  geom_line(aes(colour=Platform, linetype=Platform)) + 
  geom_ribbon(aes(group=Platform,ymin=L99,ymax=U99), alpha=.3) +
  facet_wrap( ~ z) +  
  xlab('Proportion of Hedonic Buyers (q)') + 
  ylab('Quantity Demanded') +
  # scale_shape_discrete(name="Attitude",labels=c('Utilitarian (z=0)','Hedonic (z=1)')) +
  # scale_linetype_discrete(name="Attitude",labels=c('Utilitarian (z=0)','Hedonic (z=1)')) +
  # scale_colour_discrete(name="Attitude", labels=c('Utilitarian (z=0)','Hedonic (z=1)')) +
  theme_bw()
p
ggsave('demand_by_q_z_facets_ribbon_ggplot_12.png', height=3.5,width=7.5,units='in',dpi=300)

# png('csr_demand_conf_interval_z_q_plot.png', height=6,width=8, units='in', res=300)
# df.ci <- na.omit(df.comb)
# p <- ggplot(aes(x=q,y=counts), data=df.comb) +
#   geom_point(aes(colour=Platform)) +
#   # geom_line(aes(x = q, y=mu), data=df.ci) +
#   # geom_smooth() +
#   geom_ribbon(aes(ymin=L95,ymax=U95), data=df.ci,alpha=0.9) +
#   facet_grid(z ~ Platform) +
#   xlab('Proportion of Hedonic Demand (q)') +
#   theme_bw()
# p
# dev.off()



## JAGS MODEL  

# getModelstring <- function(sig1,sig2)
# {
#   paste0("model{
#     for (i in 1:n) {
#       z[i] ~ dbern(q)
#       th1[i] <- p2*(v1 + w*sig1*z[i])
#       th2[i] <- p1*(v2 + w*sig2*z[i])
#       s[i] <- ",ifelse(sig2 > sig1,
#         "((th2[i]/th1[i])*pow(J2, rho)) / ( pow(J1, rho) + (th2[i]/th1[i])*pow(J2, rho) )",
#         "((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )"
#       ),"
#       G[i] ~ dbinom(s[i],L)
#     }
#     q ~ dbeta(h1t,h2t)
#   }")
# }


# getModelstring <- function(sig1,sig2)
# {
#   paste0("model{
#       for (i in 1:n) {
#         z[i] ~ dbern(q)
#         th1[i] <- p2*(v1 + w*sig1*z[i])
#         th2[i] <- p1*(v2 + w*sig2*z[i])
#         s[i] <- ((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )
#         G[i] ~ ",ifelse(sig2 > sig1,
#           "dbinom( 1-s[i], L)",
#           "dbinom( s[i],  L)"
#         ),"
#     }
#   q ~ dbeta(h1t,h2t)
# }")
# }

getModelstring <- function(sig1,sig2)
{
  "model{
    for (i in 1:n) {
      z[i] ~ dbern(q)
      th1[i] <- p2*(v1 + w*sig1*z[i])
      th2[i] <- p1*(v2 + w*sig2*z[i])
      s[i] <- ((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )
      G[i] ~ dbinom(s[i],L)
    }
    q ~ dbeta(h1t,h2t)
  }
  "
}

## PARAMS
q <- 0.2
N <- 500
a1 <- a2 <- 1
h1 <- 1
h2 <- 1
w <- 10
J1 <- 300
J2 <- 700
p1 <- p2 <- 10
v1 <- v2 <- 1.1
rho <- 1
sig1 <- 1
sig2 <- 0
Y <- 2000
theta <- getTheta(v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,w=ifelse(q>0,w,0),  ## ensure no signal when q=0
                  rho=rho,sig1=sig1,sig2=sig2)
L <- floor(Y/mean(theta$p1,theta$p2))

## SIMULATE DATA
df.sim <- getSimulation(q=q,N=N, theta=theta)




## z unobserved -----------------------------------------------------
data <- list(G=df.sim$G1,  L=L, 
             n=nrow(df.sim),  ## NO Z HERE
             sig1=sig1, sig2=sig2,v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,w=ifelse(q>0,w,0),    ## ensure no signal when q=0
             rho=rho,h1t=a1+h1,h2t=a2+h2)
n.chains <- 3
model <- jags.model(textConnection(getModelstring(sig1,sig2)),
                    data=data,n.adapt=3000,n.chains=n.chains )
update(model, n.iter=1000)
output <- coda.samples(model=model, variable.names=c('q'),n.iter=3000, thin=3)

fig.name <- paste0('mcmc_diagnostics_UN_q_',q,'_w_',w,'_sig1_',sig1,'_sig2_',sig2,'.png')
png(fig.name,height=6,width = 8, units = 'in', res = 250)
  par(mfrow=c(2,2),mar=c(4,3,2,2))
  mcmcplots::denplot(output,style = 'plain',auto.layout = F,main="Density of q")
  mcmcplots::traplot(output,style = 'plain',auto.layout = F,main="Trace of q")
  mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main="Thinned Running Mean of q")
  mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main="Autocorrelation of q")
dev.off()

(su <- summary(output))
(ci <- su$quantiles[c(1,5)])
(hdi <- c(hdi=unname(diff(ci))))
out <- list(
  hedonic_buyers=sum(df.sim$z),
  revenue=sum(data$p1 * data$G),
  cost='?',
  profit='?'
)
out


##  z observed -----------------------------------------------------------------
data <- list(G=df.sim$G1, z=df.sim$z, L=L, 
             n=nrow(df.sim)-1,
             sig1=sig1, sig2=sig2,v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,w=w,rho=rho,
             h1t=a1+h1,h2t=a2+h2)
n.chains <- 3
model <- jags.model(textConnection(modelstring),data=data,n.adapt=3000,n.chains=n.chains)
update(model, n.iter=1000)
output <- coda.samples(model=model, variable.names=c('q'),n.iter=6000, thin=3)

fig.name <- paste0('mcmc_diagnostics_obs_q_',q,'_w_',w,'_sig1_',sig1,'_sig2_',sig2,'.png')
  png(fig.name,height=6,width = 8, units = 'in', res = 250)
  par(mfrow=c(2,2),mar=c(4,3,2,2))
  mcmcplots::denplot(output,style = 'plain',auto.layout = F,main="Density of q")
  mcmcplots::traplot(output,style = 'plain',auto.layout = F,main="Trace of q")
  mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main="Thinned Running Mean of q")
  mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main="Autocorrelation of q")
dev.off()

# mcmcplots::denoverplot(output[[1]],output[[2]])

(su <- summary(output))
(ci <- su$quantiles[c(1,5)])
(hdi <- c(hdi=unname(diff(ci))))
out <- list(
  hedonic_buyers=sum(df.sim$z),
  revenue=sum(data$p1 * data$G),
  cost='?',
  profit='?'
)
out


#----------------------------------------------------------------

# modelstring2 <- "
# model{  
#   q ~ dbeta(h1t,h2t)
#   z ~ dbinom(q,n)
#   qhat <- z/n
#   th1 <- p2*(v1 + w*sig1*qhat)
#   th2 <- p1*(v2 + w*sig2*qhat)
#   s <- ((th1/th2)*pow(J1, rho)) / ( (th1/th2)*pow(J1, rho) + pow(J2, rho) )
#   G ~ dbinom(s,L)
# }
# "
# 
# q <- .5
# N <- 200
# a1 <- a2 <- 1
# h1 <- 1
# h2 <- 1
# df.sim <- getSimulation(q=q,N=N)
# data <- list(G=sum(df.sim$G1),  L=L*nrow(df.sim), 
#              n=nrow(df.sim), 
#              sig1=1, sig2=0, v1=1.1,v2=1.1,J1=300,J2=700,p1=10,p2=10,w=1.1,rho=.7,
#              h1t=a1+h1,h2t=a2+h2)
# 
# n.chains <- 3
# model <- jags.model(textConnection(modelstring2),data=data,n.adapt=2000,n.chains=n.chains)
# update(model, n.iter=1000)
# output <- coda.samples(model=model, variable.names=c('q'),n.iter=4000, thin=2)
# print(summary(output))
# plot(output)
# mcmcplots::autplot1(output, chain=n.chains)
# mcmcplots::denoverplot(output[[1]],output[[2]])
# mcmcplots::rmeanplot(output)
# mcmcplots::traplot(output)
# 
# su <- summary(output)
# su$quantiles[c(1,5)]
# 
# out <- list(
#   hedonic_buyers=sum(df.sim$z),
#   revenue=sum(data$p1 * data$G),
#   cost='?',
#   profit='?'
# )
# out

#------------------------- CSR BAYES GAME ----------------------------------
# 
#                   CHANGE ALL DATAFRAME TO LISTS
# 
#---------------------------------------------------------------------------
getGrowingVector <- function(N0,Tau,growth)
{
  out <- lapply(seq(0,Tau-1),function(t){
    n <- ceiling(N0*(1+growth)^t)
    return(rep(NA,n))
  })
  return(out)
}
getPsi <- function(gamma,y,p,B)
{
  return(gamma / ((y/p)*B))
}
getB <- function(s,m,b,d)
{
  return(s*m + b*(1-d))
}
getJ <- function(y,rho,gamma,c,B,f,J,dj)
{
  netMargProfit <- ((rho+1)/rho) - (gamma/(rho*c))
  return(y*netMargProfit*(B/f) + J*(1-dj))
}
getG <- function(s,L,M,seed=1111)
{
  set.seed(seed)
  s.sample <- sample(s,M,replace = F)
  return( sapply(s.sample, function(s)rbinom(n=1, size = L, s)) )
}
getQty <- function(y,p,netB,G)
{
  return(netB*(y/p) + G)
}
getPi <- function(r,w,psi,rho,c,Q,O,y)
{
  netMargProfit <-  r - ( (w+psi)/(rho*c) )
  return( netMargProfit*Q - (O/y) )
}


getModelstring <- function(sig1,sig2)
{
  paste0("model{
    for (i in 1:n) {
      z[i] ~ dbern(q)
      th1[i] <- p2*(v1 + w*sig1*z[i])
      th2[i] <- p1*(v2 + w*sig2*z[i])
      s[i] <- ",ifelse(sig2 > sig1,
        "((th2[i]/th1[i])*pow(J2, rho)) / ( pow(J1, rho) + (th2[i]/th1[i])*pow(J2, rho) )",
        "((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )"
      ),"
      G[i] ~ dbinom(s[i],L)
    }
    q ~ dbeta(h1t,h2t)
  }")
}

getQhatMcmc <- function(data, modelstring, variables=c('q'), n.chains=2, n.adapt=3000,n.iter.update=3000, n.iter.samples=3000, thin=3, seed=1111)
{
  cat('Starting MCMC . . . ')
  set.seed(1111)
  model <- jags.model(textConnection(modelstring),
                      data=data,n.adapt=n.adapt,n.chains=n.chains )
  update(model, n.iter=n.iter.update)
  output <- coda.samples(model=model, variable.names=variables,n.iter=n.iter.samples, thin=thin)
  return(output)
}


getQhatEst <- function(mcmc.output, probs, burninProportion=.2)
{
  samps <- unlist(sapply(mcmc.output,function(x){
      len <- length(x)
      burn <- ceiling(burninProportion*len)
      x[burn:len]
    }))
  return(quantile(samps, probs))
}

getQstarSig20 <- function(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
{
  Qty <- (y/p1) * B1
  psi1 <- gamma1 / Qty
  num <- psi1*v1*( p2*(v1+omega)*J1^rho  + p1*v2*J2^rho )
  denom <- omega*v2*p1*J2^rho * ( r1*rho*c1 - (w1+psi1) )
  return(num/denom)
}
getQstarSig21 <- function(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
{
  psi1 <- gamma1 / ((y/p1)*B1)
  num <- (psi1*v1*p2*J1^rho)/(v1*p2*J1^rho + v2*p1*J2^rho)
  den.a <- ((r1*rho*c1-(w1+psi1))*(v1+omega)*p2*J1^rho)/((v1+omega)*p2*J1^rho + v2*p1*J2^rho)
  den.b <- ((r1*rho*c1-w1)*v1*p2*J1^rho)/(v1*p2*J1^rho + (v2+omega)*p1*J2^rho)
  den.c <- (psi1*v1*p2*J1^rho)/(v1*p2*J1^rho + v2*p1*J2^rho)
  return(num/(den.a+den.b+den.c))
}

getSigmaStar <- function(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1, qhat)
{
  ## FIND OPPONENT CUTOFF q
  qstar2.if.sig1.0 <- getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
  qstar2.if.sig1.1 <- getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1)
  # FIND OPPONENT STRATEGY DECISION AS FUNTION OF OWN STRATEGY 
  sig2.if.sig1.0 <- ifelse(qstar2.if.sig1.0 > qhat, 1, 0)
  sig2.if.sig1.1 <- ifelse(qstar2.if.sig1.1 > qhat, 1, 0)
  ## KNOWING WHAT THE OPPONENT WILL PLAY IF WE PLAY 0,  GET PROFITABLE STRATEGY COMPARED TO qhat
  qstar1.of.sig2.if.sig1.0 <- ifelse(sig2.if.sig1.0==0,
                                   getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1),
                                   getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1))
  sig1.of.sig2.if.sig1.0 <- ifelse(qstar1.of.sig2.if.sig1.0 > qhat, 1, 0)
  ## KNOWING WHAT THE OPPONENT WILL PLAY IF WE PLAY 1,  GET PROFITABLE STRATEGY COMPARED TO qhat
  qstar1.of.sig2.if.sig1.1 <- ifelse(sig2.if.sig1.1==0,
                                   getQstarSig20(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1),
                                   getQstarSig21(omega,rho,r1,c1,w1,v1,v2,p1,p2,J1,J2,y,gamma1,B1))
  sig1.of.sig2.if.sig1.1 <- ifelse(qstar1.of.sig2.if.sig1.1 > qhat, 1, 0)
  return()  # ??????????????????????????????????????
}

#------------------------ CHECK BEHAVIOR OF QSTAR ------------------------------------
getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,4000)

lo <- 100
omegavec <- seq(.01,10,length.out = lo)
rhovec <- seq(.01,2.5,length.out = lo)
rvec <- seq(.001,.99, length.out = lo)
cvec <- seq(.01,.99, length.out = lo)
wvec <- seq(.01,.99, length.out = lo)
pvec <- seq(1,1000,length.out = lo)
jvec <- seq(10,10000,length.out = lo)
bvec <- seq(10,10000,length.out = lo)
png('qstar_sig2_0.png',height=7,width=7,units='in',res=250)
  par(mfrow=c(3,3),mar=c(4,4.5,2,1))
  x <- omegavec; plot(x,getQstarSig20(x,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(omega),ylab=expression(q^star));abline(h=c(0,1))
  x <- rhovec; plot(x,getQstarSig20(1.5,x,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(rho),ylab=expression(q^star));abline(h=c(0,1))
  x <- rvec; plot(x,getQstarSig20(1.5,1.05,x,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(r[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- cvec; plot(x,getQstarSig20(1.5,1.05,.1,x,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(c[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- wvec; plot(x,getQstarSig20(1.5,1.05,.1,.5,x,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(w[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- pvec; plot(x,getQstarSig20(1.5,1.05,.1,.5,.05,1,2,x,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(p[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- jvec; plot(x,getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,x,600,100,.1,4000),type='o',pch=16,xlab=expression(J[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- jvec; plot(x,getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,400,x,100,.1,4000),type='o',pch=16,xlab=expression(J[2]),ylab=expression(q^star));abline(h=c(0,1))
  x <- bvec; plot(x,getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,x),type='o',pch=16,xlab=expression(B[1]),ylab=expression(q^star));abline(h=c(0,1))
dev.off()

png('qstar_sig2_1.png',height=7,width=7,units='in',res=250)
  par(mfrow=c(3,3),mar=c(4,4.5,2,1))
  x <- omegavec; plot(x,getQstarSig21(x,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(omega),ylab=expression(q^star));abline(h=c(0,1))
  x <- rhovec; plot(x,getQstarSig21(1.5,x,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(rho),ylab=expression(q^star));abline(h=c(0,1))
  x <- rvec; plot(x,getQstarSig21(1.5,1.05,x,.5,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(r[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- cvec; plot(x,getQstarSig21(1.5,1.05,.1,x,.05,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(c[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- wvec; plot(x,getQstarSig21(1.5,1.05,.1,.5,x,1,2,100,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(w[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- pvec; plot(x,getQstarSig21(1.5,1.05,.1,.5,.05,1,2,x,100,400,600,100,.1,4000),type='o',pch=16,xlab=expression(p[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- jvec; plot(x,getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,x,600,100,.1,4000),type='o',pch=16,xlab=expression(J[1]),ylab=expression(q^star));abline(h=c(0,1))
  x <- jvec; plot(x,getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,400,x,100,.1,4000),type='o',pch=16,xlab=expression(J[2]),ylab=expression(q^star));abline(h=c(0,1))
  x <- bvec; plot(x,getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,x),type='o',pch=16,xlab=expression(B[1]),ylab=expression(q^star));abline(h=c(0,1))
dev.off()

png('qstar_sig2_0_sig2_1_compare_matplot1.png',height=8,width=9,units='in',res=250)
  par(mfrow=c(3,3),mar=c(4,4.5,2,1))
  x <- omegavec; df.x <- data.frame(
    sig20=getQstarSig20(x,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),
    sig21=getQstarSig21(x,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(omega),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- rhovec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,x,.1,.5,.05,1,2,100,100,400,600,100,.1,4000),
    sig21=getQstarSig21(1.5,x,.1,.5,.05,1,2,100,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(tilde(rho)),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- rvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,x,.5,.05,1,2,100,100,400,600,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,x,.5,.05,1,2,100,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(r[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- cvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,x,.05,1,2,100,100,400,600,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,.1,x,.05,1,2,100,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(c[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- wvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,.5,x,1,2,100,100,400,600,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,.1,.5,x,1,2,100,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(w[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- pvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,.5,.05,1,2,x,100,400,600,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,.1,.5,.05,1,2,x,100,400,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(p[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- jvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,x,600,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,x,600,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(J[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- jvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,400,x,100,.1,4000),
    sig21=getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,400,x,100,.1,4000)
  ); matplot(x,df.x,type='l',xlab=expression(J[2]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
  x <- bvec; df.x <- data.frame(
    sig20=getQstarSig20(1.5,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,x),
    sig21=getQstarSig21(1.5,1.05,.1,.5,.05,1,2,100,100,400,600,100,.1,x)
  ); matplot(x,df.x,type='l',xlab=expression(B[1]),ylab=expression(q^'*'));abline(h=c(0,1));legend('topright',legend=c('sigma2=0','sigma2=1'),col=1:2,lty=1:2)
dev.off()

#-----------------------------------------------

x <- list(
    v1=1
  , v2=1
  , db1=.3  # 30% buy all (y/pk) goods from current platform 1; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
  , db2=.3  # 30% buy all (y/pk) goods from current platform 2; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
  , dj1=.05
  , dj2=.05
  , c1=.5       ## seller MARGINAL cost
  , c2=.5       ## seller MARGINAL cost
  , gamma1=.05  ## seller CSR cost
  , gamma2=.05  ## seller CSR cost
  , d1=.01      ## Platform operator MARGINAL cost
  , d2=.01      ## Platform operator MARGINAL cost
  , psi1=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
  , psi2=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
  , a1=1
  , a2=1
  , r1=.1
  , r2=.1
  , w=1.5
  , rho=.7
  , growth=.01
  , Y=1000
  , ep=1e-1
  , N0=500
  , Tau=6
  , probs=c(.005,.025,.5,.975,.995)
  , learningThreshold=.05
  , n.iter=1000
  , downweight=FALSE
  , q=.7
)
x$N <- ceiling(x$N0*(1+x$growth)^(x$Tau-1))


## allocate game array
l <- list(
    M=rep(0,x$Tau)
  , qhat=list(mcmc=list(),
              est=data.frame(L99=rep(0,x$Tau),L95=rep(0,x$Tau),
                             mu=rep(0,x$Tau),
                             U95=rep(0,x$Tau),U99=rep(0,x$Tau)))
  , qstar=data.frame(qstar1=rep(0,x$Tau), qstart2=rep(0,x$Tau))
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
  , z=getGrowingVector(x$N0,x$Tau,x$growth)
  , s=getGrowingVector(x$N0,x$Tau,x$growth)
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
l$J$J1[t] <- 30
l$J$J2[t] <- 70
l$B$B1[t] <- 300
l$B$B2[t] <- 700
l$p$p1[t] <- 10
l$p$p2[t] <- 10
l$sig$sig1[t] <- 1
l$sig$sig2[t] <- 0
l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1,x$gamma2, 0)
l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, getPsi(l$gamma$gamma1[t],x$Y,l$p$p1[t],l$B$B1[t]), 0)
l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, getPsi(l$gamma$gamma2[t],x$Y,l$p$p2[t],l$B$B2[t]), 0)
l$M[t] <- 0 + x$db1*l$B$B1[t] + x$db2*l$B$B2[t]+ x$dj1*l$J$J1[t] + x$dj2*l$J$J2[t]
l$z[[t]] <- rbinom(length(l$z[[t]]), 1, x$q)  ##  CHANGING TO GROUND TRUE  q
l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])



# LIST demand share
l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                   x$v1, x$v2,
                   l$sig$sig1[t], l$sig$sig2[t],
                   l$J$J1[t], l$J$J2[t],
                   x$w, l$z[[t]],
                   x$rho, k=1)
l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
l$G$G2[[t]] <- getG( 1-l$s[[t]], l$L$L2[t], l$M[t])

## Qstar THRESHOLD
# qstar1.0 <- getQstarSig20(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t])
# qstar1.1 <- getQstarSig21(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t])
# qstar2.0 <- getQstarSig20(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t])
# qstar2.1 <- getQstarSig21(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t])
# l$qstar$qstar1[t] <- getSigmaStar(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma1,l$B$B1[t], l$qhat$est$mu[t]) ## ??????????????????
# l$qstar$qstar2[t] <- getSigmaStar(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t],x$Y,x$gamma2,l$B$B2[t], l$qhat$est$mu[t]) ## ??????????????????

## LEARN Qhat
data <- list(G=l$G$G1[t][[1]],  L=l$L$L1[t], 
             n=round(l$M[t]),  ## NO Z HERE
             sig1=l$sig$sig1[t], sig2=l$sig$sig2[t],
             J1=l$J$J1[t],J2=l$J$J2[t],p1=l$p$p1[t],p2=l$p$p2[t],
             v1=x$v1,v2=x$v2,
             w=ifelse(l$qhat$est$mu[t]>0, x$w, 0),    ## ensure no signal when q=0
             rho=x$rho,
             h1t=x$a1 + l$h$h1[t], h2t=x$a2 + l$h$h2[t])
modelstring <- getModelstring(l$sig$sig1[t], l$sig$sig2[t])
l$qhat$mcmc[[t]]  <- getQhatMcmc(data, modelstring, variables=c('q'), 
                              n.chains=3, n.adapt=300,n.iter.update=300, 
                              n.iter.samples=300, thin=3, seed=1111)
l$qhat$est[t, ] <- getQhatEst(l$qhat$mcmc[[t]], probs=x$probs, burninProportion = .2)

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
if(x$downweight) {
  ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=F, type='Ljung-Box')
  l$h$h1[t] <- x$a1 + sum(l$z[[t]])/mean(unlist(ac$statistic))
  l$h$h2[t] <- x$a2 + ( length(l$z[[t]]) - sum(l$z[[t]]) ) /mean(unlist(ac$statistic))
} else {
  l$h$h1[t] <- x$a1 + sum(l$z[[t]])
  l$h$h2[t] <- x$a2 + ( length(l$z[[t]]) - sum(l$z[[t]]) )
}


#### OUTCOME2
## Quantity
l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
## Platform Operator Profit
l$Pi$Pi1[t] <- getPi(x$r1,x$d1,l$psi$psi1[t],x$rho,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
l$Pi$Pi2[t] <- getPi(x$r2,x$d2,l$psi$psi2[t],x$rho,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)




#---------------------------------- MAIN GAME LOOP ----------------------------------------------

for (t in 2:x$Tau)
{
  cat(paste0('\nt: ',t,'\n'))
  ## STRATEGY DECISION VARIABLES
  # l$qstar$qstar1[t] <- getSigmaStar(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma1,l$B$B1[t-1], l$qhat$est$mu[t-1]) ## ??????????????????
  # l$qstar$qstar2[t] <- getSigmaStar(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma2,l$B$B2[t-1], l$qhat$est$mu[t-1]) ## ??????????????????
  l$qstar$qstar1[t] <-  ifelse(l$sig$sig2[t]==0,
                               getQstarSig20(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma1,l$B$B1[t-1]),
                               getQstarSig21(x$w,x$rho,x$r1,x$c1,x$d1,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma1,l$B$B1[t-1])
                               )
  l$qstar$qstar2[t] <- ifelse(l$sig$sig2[t]==1,
                              getQstarSig20(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma2,l$B$B2[t-1]),
                              getQstarSig21(x$w,x$rho,x$r2,x$c2,x$d2,x$v1,x$v2,l$p$p1[t-1],l$p$p2[t-1],l$J$J1[t-1],l$J$J2[t-1],x$Y,x$gamma2,l$B$B2[t-1])
                              )
  l$sig$sig1[t] <- ifelse(l$qhat$est$mu[t-1] > l$qstar$qstar1[t], 1, 0)
  l$sig$sig2[t] <- ifelse(l$qhat$est$mu[t-1] > l$qstar$qstar1[t], 1, 0)
  
  ## CSR CONTINGENT COSTS
  l$gamma$gamma1[t] <- ifelse(l$sig$sig1[t]==1, x$gamma1, 0)
  l$gamma$gamma2[t] <- ifelse(l$sig$sig2[t]==1, x$gamma2, 0)
  l$psi$psi1[t] <- ifelse(l$sig$sig1[t]==1, x$psi1, 0)
  l$psi$psi2[t] <- ifelse(l$sig$sig2[t]==1, x$psi2, 0)
  
  ## PERIOD PARAMETERS
  l$p$p1[t] <- 10
  l$p$p2[t] <- 10
  l$f$f1[t] <- 1
  l$f$f2[t] <- 1
  l$L$L1[t] <- ceiling(x$Y / l$p$p1[t])
  l$L$L2[t] <- ceiling(x$Y / l$p$p2[t])
  ## PREVIOS PERIOD DEPENDENT UPDATES
  m.temp <- x$db1*l$B$B1[t-1] + x$db2*l$B$B2[t-1]+ x$dj1*l$J$J1[t-1] + x$dj2*l$J$J2[t-1] 
  l$M[t] <- ifelse( m.temp <= length(l$z[[t]]), m.temp, length(l$z[[t]]) )
  l$J$J1[t] <- getJ(x$Y,x$rho,l$gamma$gamma1[t],x$c1,l$B$B1[t-1],l$f$f1[t],l$J$J1[t-1],x$dj1) 
  l$J$J2[t] <- getJ(x$Y,x$rho,l$gamma$gamma2[t],x$c2,l$B$B2[t-1],l$f$f2[t],l$J$J2[t-1],x$dj2)
  s1.avg <- sum(l$G$G1[[t-1]]) / (length(l$G$G1[[t-1]])*l$L$L1[t])
  s2.avg <- sum(l$G$G2[[t-1]]) / (length(l$G$G2[[t-1]])*l$L$L2[t])
  l$B$B1[t] <- getB(s1.avg, l$M[[t]], l$B$B1[t], x$db1)         
  l$B$B2[t] <- getB(s2.avg, l$M[[t]], l$B$B1[t], x$db1) 
  
  # SAMPLE MARKET ATTITUDES
  l$z[[t]] <- rbinom(length(l$z[[t]]),1, x$q)
  
  # LIST demand share
  l$s[[t]] <- share(l$p$p1[t], l$p$p2[t], 
                    x$v1, x$v2,
                    l$sig$sig1[t], l$sig$sig2[t],
                    l$J$J1[t], l$J$J2[t],
                    x$w, l$z[[t]],
                    x$rho, k=1)
  l$G$G1[[t]] <- getG(l$s[[t]], l$L$L1[t], l$M[t])
  l$G$G2[[t]] <- getG(1-l$s[[t]], l$L$L2[t], l$M[t])
  
  
  ## LEARN Qhat
  data <- list(G=l$G$G1[t][[1]],  L=l$L$L1[t], 
               n=round(l$M[t]),  ## NO Z HERE
               sig1=l$sig$sig1[t], sig2=l$sig$sig2[t],
               J1=l$J$J1[t],J2=l$J$J2[t],p1=l$p$p1[t],p2=l$p$p2[t],
               v1=x$v1,v2=x$v2,
               w=ifelse(l$qhat$est$mu[t]>0, x$w, 0),    ## ensure no signal when q=0
               rho=x$rho,
               h1t=x$a1 + l$h$h1[t-1], h2t=x$a2 + l$h$h2[t-1])
  modelstring <- getModelstring(l$sig$sig1[t], l$sig$sig2[t])
  l$qhat$mcmc[[t]] <- getQhatMcmc(data, modelstring, variable=c('q'), 
                                n.chains=3, n.adapt=1000,n.iter.update=x$n.iter, 
                                n.iter.samples=x$n.iter, thin=3, seed=1111)
  l$qhat$est[t, ] <- getQhatEst(l$qhat$mcmc[[t]], probs=x$probs, burninProportion = .2)

  ## PLOTS MCMC DIAGNOSTICS
  plotMCMCdiagnostics(l$qhat$mcmc[[t]])
  
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
  if(x$downweight) {
    ac <- autocorrTestMcmc(l$qhat$mcmc[[t]], nlags=20, pvalOnly=F, type='Ljung-Box')
    l$h$h1[t] <- l$h$h1[t-1] + sum(l$z[[t]])/mean(unlist(ac$statistic))
    l$h$h2[t] <- l$h$h2[t-1] + ( length(l$z[[t]]) - sum(l$z[[t]]) ) /mean(unlist(ac$statistic))
  } else {
    l$h$h1[t] <- l$h$h1[t-1] + sum(l$z[[t]])
    l$h$h2[t] <- l$h$h2[t-1] + ( length(l$z[[t]]) - sum(l$z[[t]]) ) 
  }

  #### OUTCOME2
  ## Quantity
  l$Q$Q1[t] <- getQty(x$Y, l$p$p1[t], l$B$B1[t]*(1-x$db1), sum(l$G$G1[[t]]))
  l$Q$Q2[t] <- getQty(x$Y, l$p$p2[t], l$B$B2[t]*(1-x$db2), sum(l$G$G2[[t]]))
  ## Platform Operator Profit
  l$Pi$Pi1[t] <- getPi(x$r1,x$d1,l$psi$psi1[t],x$rho,x$c1,l$Q$Q1[t], l$O$O1[t], x$Y)
  l$Pi$Pi2[t] <- getPi(x$r2,x$d2,l$psi$psi2[t],x$rho,x$c2,l$Q$Q2[t], l$O$O2[t], x$Y)
}
#---------------------------------------------------------------------------------

print(l$sig)
matplot(l$J,type='l')
matplot(l$B/rowSums(l$B),type='l')
matplot(l$Q/rowSums(l$Q),type='l')
matplot(l$Pi/rowSums(l$Pi),type='l')
matplot(l$qhat$est,type='l',ylim=c(0,1))

# for (t in 1:x$Tau) {
#   output <- l$qhat$mcmc[[t]]
#   n.chains <- length(l$qhat$mcmc[[t]])
# # fig.name <- paste0('mcmc_diagnostics_obs_q_',q,'_w_',w,'_sig1_',sig1,'_sig2_',sig2,'_t_',t,'.png')
#   # png(fig.name,height=6,width = 8, units = 'in', res = 250)
#   par(mfrow=c(2,2),mar=c(4,3,2,2))
#   mcmcplots::denplot(output,style = 'plain',auto.layout = F,main="Density of q")
#   mcmcplots::traplot(output,style = 'plain',auto.layout = F,main="Trace of q")
#   mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main="Thinned Running Mean of q")
#   mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main="Autocorrelation of q")
# # dev.off()
# }

plotMCMCdiagnostics <- function(output)
{
  n.chains <- length(output)
  par(mfrow=c(2,2),mar=c(4,3,2,2))
  mcmcplots::denplot(output,style = 'plain',auto.layout = F,main="Density of q")
  mcmcplots::traplot(output,style = 'plain',auto.layout = F,main="Trace of q")
  mcmcplots::rmeanplot(output,style = 'plain',auto.layout = F,main="Thinned Running Mean of q")
  mcmcplots::autplot1(output, chain=n.chains,style = 'plain', main="Autocorrelation of q")
}











