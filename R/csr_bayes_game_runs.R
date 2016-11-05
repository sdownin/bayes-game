##
#   CSR BAYES GAME TESTING RUNS
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_main.R'))


##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST X

t1.change.pd <- 3              # platform 1 adds CSR policy at period
Tau <- 5                       # number of periods

## GAME CONFIG
x <- list(t=1
  , q= .4           # focal parameter
  , epsilon = .5   # focal parameter
  , J1.0=10, J2.0=20  # secondary focal param
  , p1.0=10, p2.0=10
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
  , Y=10
  , ep=1e-1
  , Tau=Tau
  , probs=c(.005,.025,.5,.975,.995)
  , learningThreshold=.05
  , n.iter=1000
  , sig1=c(rep(0,t1.change.pd),rep(1,Tau-t1.change.pd))
  , sig2=rep(1,Tau)
  , t1.change=t1.change.pd, t2.change=0
)


##------------- RUN --------------

l <- playCsrBayesGame(x, learn=TRUE)

##---------------------------------







## OUTPUT
print(l$sig)
getCsrBayesGameSummaryPlots(x,l)
matplot(t(sapply(l$sim,function(x)x$est$q)), type='l', lty=c(3,2,1,2,3),col=c(1,1,2,1,1))
matplot(t(sapply(l$sim,function(x)x$est$epsilon)),type='l', lty=c(3,2,1,2,3),col=c(1,1,2,1,1))


## PLOT GAME SUMMARY FOR GAME WITH CSR BREAKPOINT (START CSR AT t_star)
## 
file.title <- sprintf('csr_advantage_q_%s_omeg_%s_rho_%s_w1_%s_c1_%s_downw_%s_db1_%s.png',
                      x$q,x$omega,x$rho,x$w1,x$c1,x$downweight,x$db1)
png(file.path(getwd(),'img',file.title),height=8,width=6.5,units='in',res=250)
getCsrBayesGameSummaryPlots(x,l)
dev.off()








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


