##
#   CSR BAYES GAME SIMULATIONS
#
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_main.R'))
library(ggplot2)
library(reshape2)
library(lattice); library(latticeExtra)
library(stringr)
library(foreach)
library(doParallel)

#  load(file = "l_list_7epsilons_only3finished.RData")

#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST x

## Set strategy change periods and total simulation length
t1.change.pd <- 5 # 203            # platform 1 adds CSR policy at period
t2.change.pd <- 1 # 3              # platform 2 adds CSR policy at period
Tau <- 10 # 1003                   # number of periods

## GAME CONFIG
x <- list(t=1
          , q= .5           # focal parameter
          , epsilon = 1    # focal parameter
          , params=c('q')
          , J1.0=200, J2.0=50  # secondary focal param
          #, p1.0=1, p2.0=1
          , v1= 1, v2=1
          , omega=2  ## ****
          , db1=.1, db2=.1          ## (db1)% buy all (y/pk) goods from current platform 2; (1-db1)% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
          , dj1=.1, dj2=.1
          , c1=.5, c2=.5            ## seller MARGINAL cost
          , gamma1=.2, gamma2=.2  ## seller CSR cost
          , phi1 =.4, phi2 =.4
          , w1=.02, w2=.02          ## Platform operator MARGINAL cost
          , psi1=.02, psi2=.02      ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
          , a1=1, a2=1
          , r1=.1, r2=.1
          , growth=.01
          , Y=100
          , ep=1
          , Tau=Tau
          , probs=c(.005,.025,.5,.975,.995)
          , learningThreshold=.05
          , n.iter=100
          , sig1=c(rep(0,t1.change.pd),rep(1,Tau-t1.change.pd))
          , sig2=c(rep(0,t2.change.pd),rep(1,Tau-t2.change.pd))   #rep(1,Tau)
          , t1.change=t1.change.pd, t2.change=t2.change.pd
          , cl.cutoff=0.7   # clustering between SS cutoff for 'learning' q
          , parallel = TRUE  
          , method = 'rjparallel'
          , n.cores = 4
          , learn = FALSE
          , Gdist = 'binomial'
)
#  learn methods:  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#                  'background', 'bgparallel' or 'snow'

#----------------- Multi-Parm Simulation  ----------------------------
qs <- c(0,.2,.4,.6,.8,1)   #seq(0,1,.1)
epsilons <- c(.3,1.05,1.8)  #
dbs <-  c(0.05, 0.1, 0.2) #c(.05,.5)
phis <- c(0.02, 0.1, 0.5)   # CSR cost-base price increase
l.list <- list()
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      for (l in 1:length(phis)) {
        cat(sprintf('q %s, epsilon %s, db %s, phi1 %s\n', qs[i],epsilons[j],dbs[k], phis[l]))
        x$q <- qs[i]
        x$epsilon <- epsilons[j]
        x$db1 <- x$db2 <- dbs[k]
        x$phi1 <- x$phi2 <- phis[l]
        l.list[[paste0('q',x$q,"_",'epsilon',x$epsilon,"_","db",x$db1,"_phi",x$phi1)]] <- playCsrBayesGame(x, verbose = FALSE)      
      }
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file  
# image.file <- sprintf('l_list_facet_plot_T_%s_w_%s_J1_%s_J2_%s_db_%s_q_%s.RData',x$Tau,x$omega,x$J1.0,x$J2.0,paste(dbs,collapse="-"),paste(qs,collapse="-"))
# save.image(image.file)


##--------- PLOT BUYER SHARE -------------------------
df <- getBasePlotDf(l.list, id.vars=c('q','epsilon','db','phi','period'))
db_i <- "0.2"   # 0.05, 0.1, 0.2
phi_i <- c("0.02","0.1","0.5")
df <- subset(df, subset=(db %in% db_i & phi %in% phi_i))

colourCount = length(unique(df$q))
db_num <- as.numeric(db_i)
db_lab <- ifelse(db_num < 0.1, 'Low',ifelse(db_num < 0.2, 'Moderate','High'))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
gglty <- rep(c(2,1), length(unique(df$phi))*length(unique(df$epsilon)) )

gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  scale_x_log10() +
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(phi ~ epsilon) + 
  geom_point(aes(pch=q), data=subset(df, period %% floor(max(df$period)/20) == 1)) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1),lty=gglty) +
  ggtitle(sprintf('%s Churn (%.2f)',db_lab,db_num)) + 
  ylab("Buyer Share") + xlab('Period') + theme_bw()
gg

file.name <- sprintf("buyer_base_share_PRICE_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_omega_%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau,db_i,x$omega)
ggsave(file.name, gg, height=8, width=12, units='in')





# ##----------- SELLER SHARE --------------------------
# df <- getSellerPlotDf(l.list, id.vars=c('q','epsilon','db','period'))
# ##
# gg <- ggplot(aes(x=period, y=value, colour=epsilon), data=df) + scale_x_log10() +
#   geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ q) +
#   geom_vline(xintercept=c(x$t1.change+2,x$t2.change+2)) +
#   ggtitle("Seller Share") + ylab("Seller Share") + theme_bw()
# file.name <- sprintf("seller_share_face_db_q_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
# ggsave(file.name, gg, height=8, width=16, units='in')
# ##
# gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + scale_x_log10() +
#   geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ epsilon) +
#   geom_vline(xintercept=c(x$t1.change+2,x$t2.change+2)) +
#   ggtitle("Seller Share") + ylab("Seller Share") + theme_bw()
# file.name <- sprintf("seller_share_face_db_epsilon_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
# ggsave(file.name, gg, height=8, width=16, units='in')










