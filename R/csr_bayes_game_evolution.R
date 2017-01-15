##
#   CSR BAYES GAME TESTING RUNS
#
##
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
source(file.path(getwd(),'R','csr_bayes_game_main.R'))
library(ggplot2)
library(reshape2)
library(lattice); library(latticeExtra)
library(stringr)




#  load(file = "l_list_7epsilons.RData")

#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST X

t1.change.pd <- 5            # platform 1 adds CSR policy at period
t2.change.pd <- 1              # platform 2 adds CSR policy at period
Tau <- 10                     # number of periods

## GAME CONFIG
x <- list(t=1
          , q= .5           # focal parameter
          , epsilon = 1    # focal parameter
          , params=c('q')
          , J1.0=200, J2.0=50  # secondary focal param
          , p1.0=1, p2.0=1
          , v1= 1, v2=1
          , db1=.1, db2=.1          ## (db1)% buy all (y/pk) goods from current platform 2; (1-db1)% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
          , dj1=.1, dj2=.1
          , c1=.4, c2=.4            ## seller MARGINAL cost
          , gamma1=.1, gamma2=.1  ## seller CSR cost
          , phi1=.2, phi2=.2
          , w1=.02, w2=.02          ## Platform operator MARGINAL cost
          , psi1=.01, psi2=.01      ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
          , a1=1, a2=1
          , r1=.1, r2=.1
          , omega=1
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
          , Gdist = 'binomial'
          , n.cores = 4
          , learn = FALSE
)
#  methods:  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#  'background', 'bgparallel' or 'snow'


## Multi-Parm RUN 
qs <- seq(0,1,.1)
epsilons <- c(0,.9,1,1.1,1.5)
dbs <- c(.5)
l.list <- list()
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      cat(sprintf('q %s, epsilon %s, db %s\n', qs[i],epsilons[j],dbs[k]))
      x$q <- qs[i]
      x$epsilon <- epsilons[j]
      x$db1 <- x$db2 <- dbs[k]
      l.list[[paste0('q',x$q,"_",'epsilon',x$epsilon,"_","db",x$db1)]] <- playCsrBayesGame(x, verbose = FALSE)      
    }
  }
}

file.name <- sprintf('l_list_evolution_J1_%s_J2_%s_epsilons_%s.RData', x$J1.0,x$J2.0,paste(epsilons,collapse="-"))
save.image(file.name)




##--------- BUYER SHARE -------------------------
df <- getBasePlotDf(l.list, id.vars=c('q','epsilon','db', 'period'))
##
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + scale_x_log10() +
  geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ epsilon) + 
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1)) +
  ggtitle("Buyer Base Share") + ylab("Buyer Share") + theme_bw()
gg
file.name <- sprintf("buyer_base_share_face_db_3_logx_epsilon_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
ggsave(file.name, gg, height=8, width=16, units='in')
# ##
# gg <- ggplot(aes(x=period, y=value, colour=epsilon), data=df) + scale_x_log10() +
#   geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ q) + 
#   geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1)) +
#   ggtitle("Buyer Base Share") + ylab("Buyer Share") + theme_bw()
# file.name <- sprintf("buyer_base_share_face_db_q_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
# ggsave(file.name, gg, height=8, width=16, units='in')


##----------- SELLER SHARE --------------------------
df <- getSellerPlotDf(l.list, id.vars=c('q','epsilon','db','period'))
##
gg <- ggplot(aes(x=period, y=value, colour=epsilon), data=df) + scale_x_log10() +
  geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ q) +
  geom_vline(xintercept=c(x$t1.change+2,x$t2.change+2)) +
  ggtitle("Seller Share") + ylab("Seller Share") + theme_bw()
file.name <- sprintf("seller_share_face_db_q_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
ggsave(file.name, gg, height=8, width=16, units='in')
##
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + scale_x_log10() +
  geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ epsilon) +
  geom_vline(xintercept=c(x$t1.change+2,x$t2.change+2)) +
  ggtitle("Seller Share") + ylab("Seller Share") + theme_bw()
file.name <- sprintf("seller_share_face_db_epsilon_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
ggsave(file.name, gg, height=8, width=16, units='in')










# ### IMAGE contourmap
# df <- ldply(l.list, function(l){
#   len <- nrow(l$B)
#   df <- data.frame(share=(l$B$B1[len] / sum(l$B[len,])))
#   return(df)
# }, .id="param_value")
# 
# ## Create df mat
# df <- separateParamCols(df, "param_value")
# 
# df_wide <- reshape2::dcast(df, epsilon ~ q, fun.aggregate = mean, value.var = "share")
# mat <- as.matrix(df_wide[,-1]); rownames(mat) <- df_wide[,1]
# 
# epsilon <- as.numeric(rownames(mat))
# q <- as.numeric(colnames(mat))
# image(epsilon, q, mat )
# contour(epsilon, q, mat, xlim=range(q), ylim=range(epsilon), add=T )
# 
# 
# ## Matplot of indirect network effects threshold
# df_long <- melt(mat, varnames = c('epsilon','q','db'), value.name = "share")
# df_long$epsilon <- as.factor(df_long$epsilon)
# df_long <- subset(df_long, subset=(epsilon %in% c(0,.2,.4,.6,.8,.1,1.2,1.4,1.6,1.8,2)))
# legend.title <- expression("Indirect\nNetwork\nEffexcts ("*epsilon*")")
# ggplot(aes(x=q,y=share,colour=epsilon), data=df_long) + 
#   geom_line(aes(colour=epsilon, lty=epsilon), lwd=1.1) +
#   geom_point() + theme_bw() + 
#   ggtitle("Equilibrium Demand Share") + 
#   xlab(expression("Hedonic Demand Proportion ("~q~")")) + ylab("Demand Share") + 
#   labs(colour=legend.title,lty=legend.title)
# ggsave("equilibrium_demand_share_by_q_lines_q.png", height=5,width=6,units='in')
