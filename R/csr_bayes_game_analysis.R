##
#
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_functions.R'))

#------------------ SIMULATE FROM GROUND TRUTH -----------------------------


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
runs <- 50
N <- 300
Y <- 500
qvec <- seq(0,1,.025)
omegavec <- c(.5,1,2)
rhovec <- c(.7,1,1.3)
df.comb <- data.frame()
li.comb <- list()
set.seed(111)
for (run in 1:runs) {
  df <- data.frame() 
  for(omega in omegavec) {
    for (rho in rhovec) {
       for (q in qvec) {
        df.sim <- getSimulation(q=q, N=N, Y=Y, theta = getTheta(omega=omega, rho=rho), setSeed = FALSE)
        df.sim$q <- q
        df.sim$omega <- omega
        df.sim$rho <- rho
        df <- rbind(df,df.sim)
      }
    }
  }
  df$z <- as.factor(paste0(ifelse(df$z==1,'Hedonic ','Utilitarian '),'(z=',df$z,')'))
  df.2 <- rbind(
    data.frame(Demand=as.integer(df$G1),
               z=df$z,
               q=df$q,
               Platform='1 (CSR)',
               Strength=df$omega,
               Network=df$rho),
    data.frame(Demand=as.integer(df$G2),
               z=df$z,
               q=df$q,
               Platform='2 (NO CSR)',
               Strength=df$omega,
               Network=df$rho)
  )
  df.counts <- ddply(.data = df.2, .variables = .(z,q,Platform,Strength,Network),summarise,counts=sum(Demand))
  li.comb[[run]] <- df.counts
  df.comb <- rbind(df.comb,df.counts)
}
# values=c(rgb(.8,.8,.5,.7), rgb(.1,.1,.5,.7))

# GET RIBBON 95% CREDIBLE INTERVALS
df.ci <- data.frame(z=NA, q=NA, Platform=NA, Strength=NA, Network=NA, L99=NA, mu=NA, U99=NA)
counter <- 1
for (i in levels(df.counts$z)) {
  for (j in unique(df.counts$q)) {
    for (k in levels(df.counts$Platform)) {
      for (l in unique(df.counts$Strength)) {
        for (m in unique(df.counts$Network)) {
          x <- df.comb[which(df.comb$z==i & df.comb$q==j & df.comb$Platform==k & df.comb$Strength==l & df.comb$Network==m), 'counts']
          df.ci[counter,c('z','q','Platform','Strength','Network')] <- c(i,j,k,l,m)
          df.ci[counter,c('L99','mu','U99')] <- quantile(na.omit(x), probs = c(.005, .5, .995))
          counter <- counter + 1          
        }
      }
    }
  }
}
p <- ggplot(aes(x=as.numeric(q),y=mu), data=na.omit(df.ci)) +
  # geom_smooth(aes(x=q, y=counts, colour=z), data=df.comb) +
  geom_point(aes(colour=Platform, shape=Platform)) +
  geom_line(aes(colour=Platform, linetype=Platform)) + 
  geom_ribbon(aes(group=Platform,ymin=L99,ymax=U99), alpha=.3) +
  facet_grid( Strength ~ z + Network ) +  
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


#--------------------------------------------------------------------------
## PLOT MAX PLATFORM CSR COST FOR POOLING CSR EQ AS FUNCTION OF CSR RESPONSE STRENGTH


# psibar <- x$r1*l$p$p1[t]-x$w1
# sq <- getMaxPsi(omegavec,x$rho,x,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t])
# plot(sq ~ omegavec, ylim=c(0,.021), type='o',pch=16,xlab=expression(omega),ylab=expression(psi[1]))
# abline(h=psibar, lty=2)

t <- 2
omegavec <- 2^seq(-5,5)
psibar <- x$r1*l$p$p1[t]-x$w1
rhovec <- c(0,1,2,3)
sq <- sapply(rhovec,function(r)getMaxPsi(omegavec,r,x,l$p$p1[t],l$p$p2[t],l$J$J1[t],l$J$J2[t]))
file.title <- sprintf('psi_func_omega_LnLn_pool_sep_j1_%s_j2_%s_p1_%s_p2_%s.png',l$J$J1[t],l$J$J2[t],l$p$p1[t],l$p$p2[t])
png(file.path(getwd(),'img',file.title),height=5,width=6,units='in',res=250)
  par(mar=c(4.1,4.1,1.5,1.5))
  matplot(x= omegavec, y=sq,log='xy',type='o',col='black',pch=17:20,xlab=expression(ln(omega)),
       ylab=expression(ln(psi[1])))
  legend('bottomright',title=expression(tilde(rho)),legend=rhovec,col='black',lty=1:4,pch=17:20)
  mtext(expression(bar(psi[1])),side=2,at = psibar, line=2.5, col='darkred')
  mtext(expression(Pooling:~sigma[1]*' = '*sigma[2]*' = 1'),side=3,at=8,line=-10,col='darkblue')
  mtext(expression(Separating:~sigma[1]!=sigma[2]),side=3,at=.1,line=-4,col='darkblue')
  abline(h=psibar, lty=2, col='red')
dev.off()

#------------------------------------------------------------------------
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












