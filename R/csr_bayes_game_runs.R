##
#   CSR BAYES GAME TESTING RUNS
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_main.R'))
library(ggplot2)

#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST X

t1.change.pd <- 8              # platform 1 adds CSR policy at period
Tau <- 10                       # number of periods

## GAME CONFIG
x <- list(t=1
  , q= .27           # focal parameter
  , epsilon = 1.5    # focal parameter
  , params=c('q')
  , J1.0=8, J2.0=24  # secondary focal param
  , p1.0=1, p2.0=1
  , v1= 1, v2=1
  , db1=.5, db2=.5          ## 30% buy all (y/pk) goods from current platform 2; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
  , dj1=.1, dj2=.1
  , c1=.5, c2=.5            ## seller MARGINAL cost
  , gamma1=.05, gamma2=.05  ## seller CSR cost
  , w1=.02, w2=.02          ## Platform operator MARGINAL cost
  , psi1=.02, psi2=.02      ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
  , a1=1, a2=1
  , r1=.1, r2=.1
  , omega=1
  , growth=.01
  , Y=1000
  , ep=1e-1
  , Tau=Tau
  , probs=c(.005,.025,.5,.975,.995)
  , learningThreshold=.05
  , n.iter=4000
  , sig1=c(rep(0,t1.change.pd),rep(1,Tau-t1.change.pd))
  , sig2=rep(1,Tau)
  , t1.change=t1.change.pd, t2.change=0
  , cl.cutoff=0.7   # clustering between SS cutoff for 'learning' q
  , parallel = FALSE  # problems with jags.parallel inits
  , n.cores = 4
)


## RUN 
l <- playCsrBayesGame(x, learn=TRUE)

par(mfrow=c(3,3)); for (i in 1:length(l$G$G1))hist(l$G$G1[[i]], col='lightgray')



##----------- Test cluster identification by q size -------
# qis <- rep(seq(0.005,0.4,.01),each=2)
# l.list <- list()
qis <- rep(seq(0.4,0.5,.01),each=2)
x$sig1 <- c(0,1)
x$sig2 <- c(1,1)
x$Tau <- length(x$sig1)
x$J1.0 <- 5
x$J2.0 <- 15
for (j in seq_along(qis) ) {
  i <- length(l.list) + 1
  x$q <- qis[j]
  l.list[[i]] <- playCsrBayesGame(x, learn=T)
}

q_long <- c()
bss <- c(sapply(l.list, function(l){
  sapply(seq_along(l$bss),function(i)l$bss[i])
}))
diffs <- c(sapply(l.list, function(l){
  sapply(seq_along(l$sig$sig1), function(i){
    q_long <<- c(q_long, l$q.true)
    return( ifelse(l$sig$sig1[i]==l$sig$sig2[i],0,1) )
  })
}))

df <- data.frame(q=q_long,bss=bss)
df$diffs <- factor(diffs)
#  Plot clustering accuracy
pl1 <- ggplot(aes(x=q, y=bss), data=df) + geom_hline(yintercept=0.8) + 
  geom_point(aes(colour=diffs)) +  geom_smooth(aes(colour=diffs)) + 
  ylim(0,1) + ggtitle('Buyer Type Identification Accuray by Hedonic Proportion') + 
  ylab('Between Cluster Sum of Squares') + theme_bw()
pl1
ggsave('buyer_type_id_accuracy_2.png', plot=pl1, width=6, height=5, units='in')
##---------------------------------------------------------

#------------ identify demand --------------------------
t <- 1

## assume known 2 cluster (2 states of nature)
k <- 2

cl <- getDemandClusters(l$G$G1[[t]], k=2)
plotDemandGroups(l$G$G1[[t]],cl)
Z <- getZ(x,cl,t)
h <- sum(Z)
qhat <- h / length(Z)


cl <- getDemandClusters(l, t=1, k='auto', maxcl = 5, all=F)
meanQtyVec <- cl[[k]]$cluster
for (i in 1:k)
  meanQtyVec[which(meanQtyVec==i)] <- cl[[k]]$centers[i]
plot(sapply(cl, function(x)sum(x$withinss)),type='b')

# demand clusters
plyr::count(l$s[[t]])
plyr::count(l$s[[length(l$s)]])
# hist(c(l$G$G1[[t]],l$G$G2[[t]]), breaks=20)
hist(l$G$G1[[t]], breaks=10); abline(v=cl[[k]]$centers, add=T, col='red', lwd=2)

##-------------------------------------------------------------------
# check Bayesian t-test (?)
tmp <- l$sim[[1]]$mcmc[[1]]
cor.test(tmp[,'q'], tmp[,'epsilon']) #[c('statistic','estimate','p.value','conf.int')]


## OUTPUT------------------------------------------------------------
print(l$sig)
getCsrBayesGameSummaryPlots(x,l)
matplot(t(sapply(l$sim,function(x)x$est$q)), type='l', lty=c(3,2,1,2,3),col=c(1,1,2,1,1))
matplot(t(sapply(l$sim,function(x)x$est$epsilon)),type='l', lty=c(3,2,1,2,3),col=c(1,1,2,1,1))


# ## PLOT GAME SUMMARY FOR GAME WITH CSR BREAKPOINT (START CSR AT t_star)
# ## 
# file.title <- sprintf('csr_advantage_q_%s_omeg_%s_rho_%s_w1_%s_c1_%s_downw_%s_db1_%s.png',
#                       x$q,x$omega,x$rho,x$w1,x$c1,x$downweight,x$db1)
# png(file.path(getwd(),'img',file.title),height=8,width=6.5,units='in',res=250)
# getCsrBayesGameSummaryPlots(x,l)
# dev.off()





#-------------------- SHARE HEATMAP ------------------------
##  s = f(q,epsilon) 
ep.vec <- seq(0,2,.02)
q.vec <- seq(0,1,.01)
t <- 1
share.mesh <- Vectorize(function(q,eps) {
  return(share(p1=10,p2=10,v1=1,v2=1,sig1=1,sig2=0,J1=100,J2=1000,omega=1,z=q,epsilon=eps,k=1))
})

q <- q.vec
epsilon <- ep.vec
s <- outer(q.vec, ep.vec, share.mesh)

persp.withcol(q,epsilon,s,heat.colors,20, 
              theta=-60, phi=25, expand=0.5,
              ltheta=120, ticktype='detailed', border=NA, lty=6)

image(q,epsilon,s)
k <- 11
my.cols <- brewer.pal(k, "RdYlBu")
contour(q,epsilon,s , drawlabels=T, nlevels=k, col=my.cols, add=T, lwd=2)

# per <- persp(q,epsilon,s,
#       theta=50, phi=25, expand=0.5,
#       col='lightblue', shade=0.75, ltheta=120,
#       ticktype='detailed')
#------------------------------------------------------------





# t <- seq_len(nrow(l$qhat$est))
# g <- ggplot(aes(y=mu),data=l$qhat$est) + 
#   geom_point(aes(x=t,y=mu),size=2.5)+ 
#   geom_line(aes(x=t,y=mu), lwd=1.1) + 
#   geom_ribbon(aes(x=t,ymin=L95,ymax=U95),alpha=.20) +
#   geom_ribbon(aes(x=t,ymin=L99,ymax=U99),alpha=.12) + 
#   geom_hline(yintercept=l$q.true, col='red', lty=2) +
#   ylab('q') + ylim(0,1) +
#   theme_bw()
# g

s <- 0.0054
ep <- seq(0,2,.01)
q <- getQfromEpsilonS(ep,s,sig1=0,sig2=1,J1=5,J2=100, omega=10)
plot(q, ep, xlim=c(-.01,1.01), type='l',col='red');abline(v=c(0,1))

plot(q ~ ep, ylim=c(-.01,2),type='l',col='red');abline(h=c(0,1))

#------------
N <- sum(c(2,20)*4)
h <- getQfromEpsilonS(1.5,.6,sig1=0,sig2=1,J1=2,J2=20, omega=10) * N
N - abs(h)

s.z0 <- share(p1=10,p2=10,v1=1,v2=1,sig1=0,sig2=1,J1=2,J2=20,omega=10,z=0,epsilon=1.5)
s.z1<- share(p1=10,p2=10,v1=1,v2=1,sig1=0,sig2=1,J1=2,J2=20,omega=10,z=1,epsilon=1.5)

s.delta <- s.z0 - s.z1
getQfromEpsilonS(1.5,s.delta,sig1=0,sig2=1,J1=2,J2=20, omega=10)


s <- seq(0,1,.001)
q <- getQfromEpsilonS(1.5,s,sig1=0,sig2=1,J1=5,J2=100, omega=10)
plot(q ~ s, type='l', xlim=c(0,.02)); abline(h=c(0,1));abline(v=0)
getQfromEpsilonS(1.1,mean(l$s[[1]]),sig1=0,sig2=1,J1=5,J2=100, p1=1.25,p2=1.25, omega=10)
