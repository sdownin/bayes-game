library(rjags)
library(ggplot2)
library(coda)
library(mcmcplots)


#--------------------------- FUNCTIONS -----------------------------------------

##
# 
# @argument x [list]; x$z can be n-length vector of {1,0}
# @returns [vector]  demand share for each x$z
##
share <- function(x,k=1) 
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

##
# 
# @returns list to be used as `x` argument in demand share() function
##
getTheta <- function(v1 = 1.1, v2 = 1.1, J1 = 30, J2 = 70, p1 = 100, p2 = 100, 
                  sig1 = 1, sig2 = 0, w = 1.1 , rho = .7, z = NA)
{
  return(list(v1=v1,v2=v2,J1=J1,J2=J2,p1=p1,p2=p2,sig1=sig1,sig2=sig2,w=w,rho=rho,z=z))
}

##
#
# @returns [data.frame] simulated data for evaluating CSR Bayes game via Gibbs sampling
##
getSimulation <- function(q=.5,Y=2000,N=1000)
{
  theta <- getTheta()
  L <- floor(Y/mean(theta$p1,theta$p2))  ## number of purchases each period per person i=1,...,N
  #
  set.seed(111)
  df.sim <- data.frame( z=rbinom(n=N, size = 1, q))
  df.sim$s <- share(getTheta(z=df.sim$z))
  df.sim$G1 <- sapply(df.sim$s, function(x)rbinom(n=1, size = L, x))
  df.sim$G2 <- L - df.sim$G1 
  return(df.sim)
}

#------------------ SIMULATE FROM GROUND TRUTH q=.2-----------------------------


df.sim <- getSimulation(q=.5)

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
  df.sim <- getSimulation(q=q, N=5000)
  df.sim$q <- q
  df <- rbind(df,df.sim)
}


## FACET WRAP
## GGPLOT
df$z <- as.factor(df$z)
# p1 <- ggplot(aes(x=G1), data=df) +
#   geom_density(aes(group=z, fill=z),alpha=.5, bw=.8) +
#   facet_wrap( ~ q, nrow=2 )
# p1
p2 <- ggplot(aes(x=G1), data=df) +  
  geom_histogram(aes(group=z, colour=z, fill=z),alpha=.8, bins=10) +
  facet_wrap( ~ q, nrow=2 ) + 
  theme_bw()
p2


## JAGS MODEL  
modelstring <- "
model{
  for (i in 1:n) {
    z[i] ~ dbern(q)
    th1[i] <- p2*(v1 + w*sig1*z[i])
    th2[i] <- p1*(v2 + w*sig2*z[i])
    s[i] <- ((th1[i]/th2[i])*pow(J1, rho)) / ( (th1[i]/th2[i])*pow(J1, rho) + pow(J2, rho) )
    G[i] ~ dbinom(s[i],L)
  }
  q ~ dbeta(alpha1,alpha2)
  alpha1 ~ dnorm(a1,0.5)
  alpha2 ~ dnorm(a2,0.5)
}
"

q <- .5
N <- 100
df.sim <- getSimulation(q=q,N=N)
data <- list(G=df.sim$G1, z=df.sim$z, q=df.sim$q, n=nrow(df.sim), L=L, 
             sig1=1, sig2=0, v1=1.1,v2=1.1,J1=300,J2=700,p1=10,p2=10,w=1.1,rho=.7,
             a1=1,a2=1)

n.chains <- 3
model <- jags.model(textConnection(modelstring),data=data,n.adapt=2000,n.chains=n.chains)
update(model, n.iter=1000)
output <- coda.samples(model=model, variable.names=c('q','alpha1','alpha2'),n.iter=4000, thin=2)
print(summary(output))
plot(output)
mcmcplots::autplot1(output, chain=n.chains)
mcmcplots::denoverplot(output[[1]],output[[2]])
mcmcplots::rmeanplot(output)
mcmcplots::traplot(output)

su <- summary(output)
su$quantiles['q',c(1,5)]

out <- list(
  hedonic_buyers=sum(df.sim$z),
  revenue=sum(data$p1 * data$G),
  cost='?',
  profit='?'
)
out


#------------------------- CSR BAYES GAME ----------------------------------
li <- list(
    a1=1.1
  , a2=1.1
  , w=1.1
  , rho=.7
  , growth=.01
  , Y=1000
)

N0 <- 1000
Tau <- 48
N <- ceiling(1000*(1+li$growth)^Tau)
df <- data.frame(
  p1=rep(10,N), 
  p2=rep(10,N)
)
df$J1 <- 0
df$J2 <- 0
df$B1 <- 0
df$B2 <- 0
df$M <- NA
df$z <- NA
df$g <- NA
df$G <- NA
df$sig1 <- NA
df$sig2 <- NA






























