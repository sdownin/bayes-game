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
t2.change.pd <- 2                   # platform 2 adds CSR policy at period
t1.change.pd <- 100 + t2.change.pd  # platform 1 adds CSR policy at period
Tau <- 1000 + t2.change.pd          # number of periods

## GAME CONFIG
x <- list(t=1
          , q= .5              ## hedonic probability by nature (proportion)
          , epsilon = 1        ## indirect network effects
          , params=c('q')      ## parmeters to be 'learned' by gibbs sampling
          , J1.0=200, J2.0=50  ## Initial sellers on platforms
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

l.list <- list()

qs <- c(0,.2,.4,.6,.8,1)   #seq(0,1,.1)
epsilons <- c(.5,1.05,1.6)  #
dbs <-  c(0.05, 0.1, 0.2) #c(.05,.5)
phis <- c(0.02, 0.1, 0.4)   # CSR cost-base price increase

for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      for (l in 1:length(phis)) {
        params <- list(q=qs[i], 
                    epsilon=epsilons[j], 
                    db=dbs[k], 
                    phi=phis[l])
        print(sapply(names(params),function(name)sprintf("%s",params[name])))
        x$q <- qs[i]
        x$epsilon <- epsilons[j]
        x$db1 <- x$db2 <- dbs[k]
        x$phi1 <- x$phi2 <- phis[l]
        index <- paste0("q",x$q,"_epsilon",x$epsilon,"_db",x$db1,"_phi",x$phi1)
        l.list[[index]] <- playCsrBayesGame(x, verbose = FALSE)  
        l.list[[index]]$params <- params   
      }
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file
# ## load saved binary file as follows:
# ## load('filename.RData')
image.file <- sprintf('l_list_facet_plot_T_%s_w_%s_J1_%s_J2_%s_db_%s_q_%s.RData',x$Tau,x$omega,x$J1.0,x$J2.0,paste(dbs,collapse="-"),paste(qs,collapse="-"))
save.image(image.file)



##--------- PLOT BUYER SHARE - SET ONE VALUE OF CHURN -------------------------

## subset data to plot
df <- getBasePlotDf(l.list, id.vars=c('q','epsilon','db','phi','period'))
db_i <- "0.2"   # 0.05, 0.1, 0.2
phi_i <- c("0.02","0.1","0.4")
df <- subset(df, subset=(db %in% db_i & phi %in% phi_i))

## set numerics to factors for plotting
for(var in names(l.list[[1]]$params))
  df[,var] <- as.factor(df[,var])

## subset dataframe for geom_point characters
if (length(unique(df$period))>40) {
  df.point <- subset(df, period %% floor(max(df$period)/20) == 1)
} else {
  df.point <- df
}

## prepare plotting arguments
colourCount = length(unique(df$q))
db_num <- as.numeric(db_i)
db_lab <- ifelse(db_num < 0.1, 'Low',ifelse(db_num < 0.2, 'Moderate','High'))
gglty <- rep(c(2,1), length(unique(df$phi))*length(unique(df$epsilon)) )

# labeller = list(epsilon=sapply(epsilons, function(x)sprintf("epsilon=%s",x)),
#                 phi=sapply(phis, function(x)sprintf("phi=%s",x)))
# names(labeller$epsilon) <- epsilons
# names(labeller$phi) <- phis
labels.epsilon <- sapply(epsilons, function(x)sprintf("Network Effect = %s",x))
labels.phi <- sapply(phis, function(x)sprintf("CSR Price Premium = %s",x))
names(labels.epsilon) <- epsilons
names(labels.phi) <- phis
labels <- c(labels.epsilon, labels.phi)

# legend guide
legend.guide <- guide_legend(title="Hedonic\nProportion\n(q)")

## GGPLOT object
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  scale_x_log10() +
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(phi ~ epsilon, labeller = as_labeller(labels)) + 
  geom_point(aes(pch=q), data=df.point) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1),lty=gglty) +
  ggtitle(sprintf('%s Churn (%.2f)',db_lab,db_num)) + 
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) + 
  guides(color=legend.guide, group=legend.guide,lty=legend.guide,pch=legend.guide) +
  theme_bw() 
gg  ## display plot

## save plot
file.name <- sprintf("buyer_base_share_PRICE_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_omega_%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau,db_i,x$omega)
ggsave(file.name, gg, height=8, width=12, units='in')

##----------- No title figure for article ----------------------------
## GGPLOT object
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  scale_x_log10() +
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(phi ~ epsilon, labeller = as_labeller(labels)) + 
  geom_point(aes(pch=q), data=df.point) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1),lty=gglty) +
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) +
  guides(color=legend.guide, group=legend.guide,lty=legend.guide,pch=legend.guide) +
  theme_bw() 
gg  ## display plot

## save plot
file.name <- sprintf("buyer_base_share_PRICE_6p5-10_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_omega_%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau,db_i,x$omega)
ggsave(file.name, gg, height=6.5, width=10, units='in')
##------------------------------------------------------------



##--------- PLOT BUYER SHARE - SET ONE VALUE OF CSR PRICE PREMIUM -------------------------

## subset data to plot
df <- getBasePlotDf(l.list, id.vars=c('q','epsilon','db','phi','period'))
db_i <- c("0.05", "0.1", "0.2")
phi_i <- "0.02" # c("0.02","0.1","0.5")
df <- subset(df, subset=(db %in% db_i & phi %in% phi_i))

## set numerics to factors for plotting
for(var in names(l.list[[1]]$params))
  df[,var] <- as.factor(df[,var])

## subset dataframe for geom_point characters
if (length(unique(df$period))>40) {
  df.point <- subset(df, period %% floor(max(df$period)/20) == 1)
} else {
  df.point <- df
}

## prepare plotting arguments
colourCount = length(unique(df$q))
phi_num <- as.numeric(phi_i)
phi_lab <- ifelse(phi_num < 0.1, 'Low',ifelse(phi_num < 0.2, 'Moderate','High'))
gglty <- rep(c(2,1), length(unique(df$db))*length(unique(df$epsilon)) )

## GGPLOT object
gg <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  scale_x_log10() +
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(db ~ epsilon) + 
  geom_point(aes(pch=q), data=df.point) + 
  scale_color_manual(values=colorRamps::matlab.like(colourCount)) +
  geom_vline(xintercept=c(l.list[[1]]$t1.change+1,l.list[[1]]$t2.change+1),lty=gglty) +
  ggtitle(sprintf('%s Price Premium (%.2f)',phi_lab,phi_num)) + 
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) + theme_bw()
gg  ## display plot

## save plot
file.name <- sprintf("buyer_base_share_CHURN_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_phi_%s_omega_%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau,phi_i,x$omega)
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










