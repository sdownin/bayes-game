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

#  load(file = "l_list_evolution_equilibrium.RData")
#-------------------------------------------------------------------------------------
separateParamCols <- function(df, paramCol="param_value") {
  x <- as.character(df[ ,paramCol])
  param.df <- ldply(str_split(str_replace_all(x, "[a-zA-Z]", ""), "_"))
  names(param.df) <- unlist(str_split(str_replace_all(x[1], "[\\d.]", ""), "_"))
  df <- cbind(df[,which(!(names(df) %in% paramCol))], param.df)
  return(df)
}
getBasePlotDf <- function(l.list, id.vars = c('q','epsilon','db','period')) {
  base <- ldply(l.list, function(l) {
    l$B$period <- as.numeric(seq_len(nrow(l$B)))
    l$B$B1share <- l$B$B1 / (l$B$B1 + l$B$B2)
    return(l$B)
  }, .id = "param_value")
  df <- separateParamCols(base, "param_value")
  base_long <- melt(df, id.vars = id.vars)
  ##
  df <- subset(base_long, subset = ( !(variable %in% c('B1','B2'))))
  return(df)
}
getSellerPlotDf <- function(l.list, id.vars = c('q','epsilon','db','period')) {
  sell <- ldply(l.list, function(l) {
    l$J$period <- as.numeric(seq_len(nrow(l$J)))
    l$J$J1share <- l$J$J1 / (l$J$J1 + l$J$J2)
    return(l$J)
  }, .id = "param_value")
  df <- separateParamCols(sell, "param_value")
  sell_long <- melt(df, id.vars = id.vars)
  ##
  df <- subset(sell_long, subset = ( !(variable %in% c('J1','J2'))))
  return(df)
}
#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST X

Tau <- 503                       # number of periods
t1.change.pd <- Tau            # platform 1 adds CSR policy at period
t2.change.pd <- 3              # platform 2 adds CSR policy at period

## GAME CONFIG
x <- list(t=1
          , q= .5           # focal parameter
          , epsilon = 1    # focal parameter
          , params=c('q')
          , J1.0=200, J2.0=200  # secondary focal param
          , p1.0=1, p2.0=1
          , v1= 1, v2=1
          , db1=.5, db2=.5          ## (db1)% buy all (y/pk) goods from current platform 2; (1-db1)% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
          , dj1=.1, dj2=.1
          , c1=.5, c2=.5            ## seller MARGINAL cost
          , gamma1=.05, gamma2=.05  ## seller CSR cost
          , w1=.02, w2=.02          ## Platform operator MARGINAL cost
          , psi1=.02, psi2=.02      ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
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
          , n.cores = 4
          , learn = FALSE
)
#  methods:  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#  'background', 'bgparallel' or 'snow'


## Multi-Parm RUN 
qs <- seq(0,1,.05)
epsilons <- seq(0,2,.05)
l.list <- list()
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
      cat(sprintf('q %s, epsilon %s\n', qs[i],epsilons[j]))
      x$q <- qs[i]
      x$epsilon <- epsilons[j]
      l.list[[paste0('q',x$q,"_",'epsilon',x$epsilon)]] <- playCsrBayesGame(x, verbose = FALSE)      
  }
}

save.image("l_list_evolution_equilibrium_1.RData")


df <- ldply(l.list, function(l){
  len <- nrow(l$B)
  df <- data.frame(share=(l$B$B2[len] / sum(l$B[len,])))
  return(df)
}, .id="param_value")
df <- separateParamCols(df, paramCol = "param_value")
df_wide <- reshape2::dcast(df, epsilon ~ q, fun.aggregate = mean, value.var = "share")
mat <- as.matrix(df_wide[,-1]); rownames(mat) <- df_wide[,1]

png("demand_share_equilibrium_J1_50_J2_200_omega_1_2.png",height=6, width=6, units='in', res=250)
epsilon <- as.numeric(rownames(mat))
q <- as.numeric(colnames(mat))
image(q, epsilon, t(mat) , main="Equilibrium Demand Share", ylab=expression("Indirect Network Effects ("~epsilon~")"), xlab=expression("Hedonic Demand Proportion ("~q~")"))
contour(q, epsilon, t(mat), xlim=range(q), ylim=range(epsilon), add=T )
dev.off()



## Matplot of indirect network effects threshold
df_long <- melt(mat, varnames = c('epsilon','q'), value.name = "share")
df_long$q <- factor(df_long$q)
df_long <- subset(df_long, subset=(q %in% c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
legend.title <- expression("Hedonic\nDemand\nProportion ("*q*")")
ggplot(aes(x=epsilon,y=share,colour=q), data=df_long) + 
  geom_line(aes(colour=q, lty=q), lwd=1.1) +
  geom_point() + theme_bw() + 
  ggtitle("Equilibrium Demand Share") + 
  xlab(expression("Indirect Network Effects ("~epsilon~")")) + ylab("Demand Share") +
  labs(colour=legend.title,lty=legend.title)
ggsave("equilibrium_demand_share_by_epsilon_lines_q.png", height=5,width=6,units='in')


df_long <- melt(mat, varnames = c('epsilon','q'), value.name = "share")
df_long$epsilon <- as.factor(df_long$epsilon)
df_long <- subset(df_long, subset=(epsilon %in% c(0,.2,.4,.6,.8,.1,1.2,1.4,1.6,1.8,2)))
legend.title <- expression("Indirect\nNetwork\nEffexcts ("*epsilon*")")
ggplot(aes(x=q,y=share,colour=epsilon), data=df_long) + 
  geom_line(aes(colour=epsilon, lty=epsilon), lwd=1.1) +
  geom_point() + theme_bw() + 
  ggtitle("Equilibrium Demand Share") + 
  xlab(expression("Hedonic Demand Proportion ("~q~")")) + ylab("Demand Share") + 
  labs(colour=legend.title,lty=legend.title)
ggsave("equilibrium_demand_share_by_q_lines_q.png", height=5,width=6,units='in')




##--------- BUYER SHARE -------------------------
df <- getBasePlotDf(l.list, id.vars=c('q','epsilon','db','period'))
##
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + scale_x_log10() +
  geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ epsilon) + 
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1)) +
  ggtitle("Buyer Base Share") + ylab("Buyer Share") + theme_bw()
file.name <- sprintf("buyer_base_share_face_db_epsilon_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
ggsave(file.name, gg, height=8, width=16, units='in')
##
gg <- ggplot(aes(x=period, y=value, colour=epsilon), data=df) + scale_x_log10() +
  geom_point() +  geom_line(lwd=1.1) + facet_grid(db ~ q) + 
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1)) +
  ggtitle("Buyer Base Share") + ylab("Buyer Share") + theme_bw()
file.name <- sprintf("buyer_base_share_face_db_q_J1_%s_J2_%s_t1_%s_t2_%s_T%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau)
ggsave(file.name, gg, height=8, width=16, units='in')


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










