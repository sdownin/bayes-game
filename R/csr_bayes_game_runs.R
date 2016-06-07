##
#   CSR BAYES GAME TESTING RUNS
#
#
#
#   v1=1
# , v2=1
# , db1=.3  # 30% buy all (y/pk) goods from current platform 1; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
# , db2=.3  # 30% buy all (y/pk) goods from current platform 2; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
# , dj1=.05
# , dj2=.05
# , c1=.5       ## seller MARGINAL cost
# , c2=.5       ## seller MARGINAL cost
# , gamma1=.05  ## seller CSR cost
# , gamma2=.05  ## seller CSR cost
# , w1=.01      ## Platform operator MARGINAL cost
# , w2=.01      ## Platform operator MARGINAL cost
# , psi1=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
# , psi2=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
# , a1=1
# , a2=1
# , r1=.1
# , r2=.1
# , omega=1.5
# , rho=1.05
# , growth=.01
# , Y=1000
# , ep=1e-1
# , N0=500
# , Tau=10
# , probs=c(.005,.025,.5,.975,.995)
# , learningThreshold=.05
# , downweight=TRUE
# , q=.3
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_functions.R'))


##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST X

x <- list(
  v1= 1
  , v2=1
  , db1=.1  # 30% buy all (y/pk) goods from current platform 1; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
  , db2=.1  # 30% buy all (y/pk) goods from current platform 2; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
  , dj1=.02
  , dj2=.02
  , c1=.5       ## seller MARGINAL cost
  , c2=.5       ## seller MARGINAL cost
  , gamma1=.1  ## seller CSR cost
  , gamma2=.1  ## seller CSR cost
  , w1=.01      ## Platform operator MARGINAL cost
  , w2=.01      ## Platform operator MARGINAL cost
  , psi1=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
  , psi2=.01    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
  , a1=1
  , a2=1
  , r1=.1
  , r2=.1
  , omega=2
  , rho=1.2
  , growth=.001
  , Y=1000
  , ep=1e-1
  , N0=500
  , Tau=60
  , probs=c(.005,.025,.5,.975,.995)
  , learningThreshold=.05
  , n.iter=1000
  , downweight=TRUE
  , q=.5
  , sig1.fixed=NA
  , sig2.fixed=NA
  , t1.change=NA
  , t2.change=NA
)

## MAIN GAME CALL
## SET STRATEGY
x$t1.change <- 10
x$t2.change <- 25
x$sig1.fixed <- c(rep(0,x$t1.change),rep(1,x$Tau-x$t1.change))
x$sig2.fixed <- c(rep(0,x$t2.change),rep(1,x$Tau-x$t2.change))
## RUN
l <- playCsrBayesGame(x, learn=FALSE)
## OUTPUT
print(l$sig)
getCsrBayesGameSummaryPlots(x,l)


## PLOT GAME SUMMARY FOR GAME WITH CSR BREAKPOINT (START CSR AT t_star)
## 
file.title <- sprintf('csr_advantage_seller_pd_delay_q_%s_omeg_%s_rho_%s_w1_%s_c1_%s_downw_%s_db1_%s.png',
                      x$q,x$omega,x$rho,x$w1,x$c1,x$downweight,x$db1)
png(file.path(getwd(),'img',file.title),height=8,width=6.5,units='in',res=250)
  csrBayesGameSummaryPlots(x,l)
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





