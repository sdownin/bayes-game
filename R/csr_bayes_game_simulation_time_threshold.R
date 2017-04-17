##
#   CSR BAYES GAME SIMULATIONS
#
#
##

## SET YOUR WORKING DIRECTORY HERE, like the format below:
setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')
# setwd()

source(file.path(getwd(),'R','csr_bayes_game_main.R'))
library(ggplot2)
library(reshape2)
library(lattice); library(latticeExtra)
library(stringr)
library(foreach)
library(doParallel)

## NOTE: if missing package, run:
## install.packages(c('packageName1','packageName2')) ## for N packages


#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST x

## Set strategy change periods and total simulation length
t1.change.pd <- 100            # platform 1 adds CSR policy at period
t2.change.pd <- 2            # platform 2 adds CSR policy at period
Tau <- 1400                   # number of periods

## GAME CONFIG
x <- list(t=1
          , q= .5              ## hedonic probability by nature (proportion)
          , epsilon = 1        ## indirect network effects
          , params=c('q')      ## parmeters to be 'learned' by gibbs sampling
          , J1.0=100000, J2.0=10000  ## Initial sellers on platforms
          , B1.0=40000000, B2.0=2100000
          #, p1.0=1, p2.0=1    ## price  *MOVED* to be set dynamically as function of CSR policy
          , v1= 1, v2=1        ## utilitarian value
          , omega=2  ##        ## hedonic value
          , db1=.1, db2=.1     ## buyer churn:  (db1)% buy all (quantity = y/pk) goods from current platform 2; (1-db1)% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
          , dj1=.1, dj2=.1     ## seller churn
          , c1=.5, c2=.5       ## seller MARGINAL cost
          , gamma1=.2, gamma2=.2  ## seller CSR cost 
          , phi1 =.4, phi2 =.4    ## platfom CSR cost
          , w1=.02, w2=.02        ## Platform operator MARGINAL cost
          # , psi1=.02, psi2=.02  ## Platform operator CSR cost   *MOVED* --> function of (gamma, B, y, p1)
          , a1=1, a2=1        ## time constant
          , r1=.1, r2=.1      ## platform transaction fee (percent of price)
          , growth=.01        ## ** currently not using 
          , Y=100             ## budget
          , ep=1              ## downweighted below
          , Tau=Tau           ## number of periods in simulation
          , probs=c(.005,.025,.5,.975,.995)  ## quantiles of estimated parameter to keep
          , learningThreshold=.05            ## ** don't need
          , n.iter=100                       ## ** don't need
          , t1.change=t1.change.pd, t2.change=t2.change.pd    ## ** don't need
          , cl.cutoff=0.7         ## ** don't need # clustering between SS cutoff for 'learning' q
          , parallel = TRUE       ## ** don't need
          , method = 'rjparallel' ## ** don't need
          , n.cores = 4           ## ** don't need
          , learn = FALSE         ## ** don't need
          , Gdist = 'binomial'    ## ** don't need
)
#  learn methods:  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#                  'background', 'bgparallel' or 'snow'

#---------------------------------------------------------------------
#----------------- Multi-Param Simulation  ---------------------------
#---------------------------------------------------------------------

l.list <- list()   #testing:#    i = j = k = l = m = 1

qs <- c(0, .2, .4, .6, .8, 1)  # seq(0,1,.1) ##c(0,.02,.05,.1,.2,.3,.4) # c(0,.2,.4,.6,.8,1)  #seq(0,1,.1)
epsilons <- c(.3, 1.05, 1.8)  #
dbs <-  0.1 # c(0.05, 0.1, 0.2) #c(.05,.5)
phis <- 0.1  # c(0.02, 0.1, 0.5)   # CSR cost-base price increase
t1.changes <- c(30, 90, 270) # c(20, 80, 320) 

for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      for (l in 1:length(phis)) {
        for (m in 1:length(t1.changes)) {
          params <- list( q=qs[i], 
                          epsilon=epsilons[j], 
                          db=dbs[k], 
                          phi=phis[l], 
                          t1.change=t1.changes[m])
          print(sapply(names(params),function(name)sprintf("%s",params[name])))
          x$q <- qs[i]
          x$epsilon <- epsilons[j]
          x$db1 <- x$db2 <- dbs[k]
          x$phi1 <- x$phi2 <- phis[l]
          x$t1.change <- t1.changes[m]
          index <- paste0("q",x$q,"_epsilon",x$epsilon,"_db",x$db1,"_phi",x$phi1,"_t1.change",x$t1.change)
          l.list[[index]] <- playCsrBayesGame(x, verbose = F)  
          l.list[[index]]$params <- params
        }
      }
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file
# ## load saved binary file as follows:
# ## load('filename.RData')
image.file <- sprintf('_uber_l_list_facet_plot_PERIOD_T_%s_w_%s_J1_%s_J2_%s_db_%s_q_%s_t1_%s.RData',x$Tau,x$omega,x$J1.0,x$J2.0,paste(dbs,collapse="-"),paste(qs,collapse="-"),paste(t1.changes,collapse="-"))
save(l.list, file=image.file)



##--------- PLOT BUYER SHARE -------------------------

## choose params to display
db_i <- "0.1"   # "0.05", "0.1", "0.2"
phi_i <- "0.1" # c("0.02","0.1","0.5")
q_i <- c('0', '0.2', '0.4', '0.6', '0.8', '1')  # .4
t1.change_i <- t1.changes#c(20, 60, 180)


## subset data by chosen params
df <- getBasePlotDf(l.list, id.vars=c(names(l.list[[1]]$params), 'period') )
df <- subset(df, subset=(db %in% db_i & phi %in% phi_i & t1.change %in% t1.change_i))
df <- subset(df, subset=(as.character(q) %in% q_i))  ## couldn't filter 0.6 as float ?? need to use characters

## set numerics to factors for plotting
periods <- c('Early','Mid','Late')
df$t1.change <- sapply(seq_len(nrow(df)), function(i){
  if (df$t1.change[i]==t1.changes[1]) return( sprintf('1. Early Response = %s',df$t1.change[i]) )
  if (df$t1.change[i]==t1.changes[2]) return( sprintf('2. Mid Response = %s',df$t1.change[i]) )
  if (df$t1.change[i]==t1.changes[3]) return( sprintf('3. Late Response = %s',df$t1.change[i]) )
});
for(var in names(l.list[[1]]$params))
  df[,var] <- factor(df[,var], levels = as.character(sort(unique(df[,var]))))

## subset dataframe for geom_point characters
if (length(unique(df$period))>40) {
  df.point <- subset(df, period %% floor(max(df$period)/20) == 1)
} else {
  df.point <- df
}

## prepare plotting arguments
nPeriods <- length(unique(df$period))
ncols <- length(epsilons)
xintercepts2 <- l.list[[1]]$t2.change+1
xintercepts1 <- data.frame(z = rep(t1.change_i, each=ncols), 
                           t1.change = rep(levels(df$t1.change), each=ncols), 
                           epsilon = rep(epsilons, ncols) )
colourCount <- length(unique(df$q))
db_num <- as.numeric(db_i)
db_lab <- ifelse(db_num < 0.1, 'Low',ifelse(db_num < 0.2, 'Moderate','High'))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# make labeller
neLevels <- c('Low','Moderate','High')
# periods <- c('Early','Mid','Late')
labels.epsilon = sapply(seq_along(epsilons), function(x)sprintf("%s Network Effect = %3s",neLevels[x],epsilons[x]))
labels.t1.change = sapply(seq_along(levels(df$t1.change)), function(x)levels(df$t1.change)[x])  # sapply(seq_along(t1.change_i), function(x)sprintf("%s. %s Response = %s",x,periods[x],t1.change_i[x]))
names(labels.epsilon) <- epsilons
names(labels.t1.change) <- levels(df$t1.change)
labels <- c(labels.epsilon, labels.t1.change)

# legend guide
legend.guide <- guide_legend(title="Hedonic\nProportion\n(q)")

##----------- No title figure for article ----------------------------
## GGPLOT object
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(t1.change ~ epsilon, labeller = as_labeller(labels)) + 
  geom_point(aes(pch=q), data=df.point) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=xintercepts2, lty=1) +
  geom_vline(aes(xintercept=z), xintercepts1, lty=2) +
  scale_x_log10() + ylim(0,1) +
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) + 
  guides(color=legend.guide, group=legend.guide,lty=legend.guide,pch=legend.guide) +
  theme_bw() 
gg  ## display plot

## save plot
file.name <- sprintf("_uber_buyer_base_share_PERIOD_6p5-10_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_phi_%s_omega_%s.png", x$J1.0, x$J2.0,paste(x$t1.change,collapse = "-"),x$t2.change,x$Tau,db_i,phi_i,x$omega)
ggsave(file.name, gg, height=6.5, width=10, units='in')
##------------------------------------------------------------


##-----------------with title--------------------------------
## GGPLOT object
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(t1.change ~ epsilon, labeller = as_labeller(labels)) + 
  geom_point(aes(pch=q), data=df.point) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=xintercepts2, lty=1) +
  geom_vline(aes(xintercept=z), xintercepts1, lty=2) +
  scale_x_log10() + ylim(0,1) +
  ggtitle(sprintf('%s Churn (%.2f)',db_lab,db_num)) + 
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) + 
  guides(color=legend.guide, group=legend.guide,lty=legend.guide,pch=legend.guide) +
  theme_bw() 
gg  ## display plot

## save plot
file.name <- sprintf("buyer_base_share_PERIOD_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_phi_%s_omega_%s.png", x$J1.0, x$J2.0,paste(x$t1.change,collapse = "-"),x$t2.change,x$Tau,db_i,phi_i,x$omega)
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
