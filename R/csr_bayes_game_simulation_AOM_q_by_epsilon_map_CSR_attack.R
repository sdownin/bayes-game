##
#   CSR BAYES GAME SIMULATIONS
#
#   CSR ATTACK
#
##

## SET YOUR WORKING DIRECTORY HERE, like the format below:
setwd('C:\\Users\\T430\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')
# setwd()

source(file.path(getwd(),'R','csr_bayes_game_main.R'))
library(ggplot2)
library(reshape2)
library(lattice); library(latticeExtra)
library(stringr)
library(colorRamps)

#--------------------------------------------------------------------------------------
##  RUN MAIN GAME SIMULATION
##  USING GAME SETUP LIST x



## GAME CONFIG
x <- list(Tau=2200           ## number of periods in simulation
          ##-------ATTACK-----------
          , gamma2 =    .2      # CSR action (0 < gamma2 < (1-c2)) #where c is marginal cost
          , discount2 = 0    # PRICE  (0 < discount2 <= epsilon)
          , t2.change = 2    # platform 2 attacks:  PRICE policy at period t2
          ##----------------------------
          , J1.0=200, J2.0=50 #J1.0=100000, J2.0=10000 ## Initial sellers on platforms
          , B1.0=800, B2.0=200  #B1.0=40000000, B2.0=2100000
          , v1= 1, v2=1        ## utilitarian value
          , omega=80  ##        ## hedonic value
          , dj1=.1, dj2=.1     ## seller churn
          , c1=.5, c2=.5       ## seller MARGINAL cost
          , b1=0, b2=0         ## price discount     ***** PRICE ACTION *****
          , gamma1=.2, gamma2=.2  ## seller CSR cost ****** CSR ACTION ******
          , phi1 =.4, phi2 =.4    ## platfom CSR cost
          , w1=.02, w2=.02        ## Platform operator MARGINAL cost
          , a1=1, a2=1        ## time constant
          , r1=.1, r2=.1      ## platform transaction fee (percent of price)
          , growth=.01        ## ** currently not using 
          , Y=100             ## budget
          , ep=1              ## downweighted below
          , t=1               ## start at t=1
          ## SET RANGES OF PARAM VALUES FOR PLTOTING BELOW
          # , q= .5              ## hedonic probability by nature (proportion)
          # , epsilon = 1        ## indirect network effects
          #, p1.0=1, p2.0=1    ## price  *MOVED* to be set dynamically as function of CSR policy
          # , db1=.1, db2=.1     ## buyer churn:  (db1)% buy all (quantity = y/pk) goods from current platform 2; (1-db1)% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
          # , psi1=.02, psi2=.02  ## Platform operator CSR cost   *MOVED* --> function of (gamma, B, y, p1)
          ## Deprecated in current version 
          , params=c('q')      ## parmeters to be 'learned' by gibbs sampling
          , probs=c(.005,.025,.5,.975,.995)  ## ** don't need ##  quantiles of estimated parameter to keep
          , learningThreshold=.05            ## ** don't need
          , n.iter=100                       ## ** don't need
          , cl.cutoff=0.7         ## ** don't need # clustering between SS cutoff for 'learning' q
          , parallel = TRUE       ## ** don't need
          , method = 'rjparallel' ## ** don't need
          , n.cores = 4           ## ** don't need
          , learn = FALSE         ## ** don't need
          , Gdist = 'binomial'    ## ** don't need
)
#  learn methods:  'rjags', 'simple', 'interruptible', 'parallel', 'rjparallel', 
#                  'background', 'bgparallel' or 'snow'

##-----------------------------------------
## SIMULATION SETTSIN FOR ALL SIMULATIONS
##-----------------------------------------
n   <- 18  ## granularity (length of contingency vectors to simulate)
t2r <- 100  ## response timing (if the response action is CSR or PRICE)

#=====================================================================
#----------------- 1. NO RESPONSE SIMULATION ----------------------
#---------------------------------------------------------------------

##----RESPONSE ACTION (platform 2) ------------
gamma1s    <- 0 # CSR action (0 < gamma1 < (1-c1)) #where c is marginal cost
discount1s <- 0 # PRICE  (0 < discount1 <= epsilon)
t1.changes <- x$Tau  #  c(30, 90, 270) # c(20, 80, 320) 
##-----FIXED PARAMS-----------------------------
dbs <-  0.1 # sellers churn # c(0.05, 0.1, 0.2) #c(.05,.5)
##----CONTINGENCIES TO TEST---------------------
qs <- round(seq(0,1, length.out = n), 2)  # seq(0,1,.1) ##c(0,.02,.05,.1,.2,.3,.4) # c(0,.2,.4,.6,.8,1)  #seq(0,1,.1)
epsilons <- round(seq(0.2,2, length.out = n), 2)  #
##----init data objects-------------------------
l.list <- list()  
z <- matrix(NA, nrow = length(epsilons), ncol=length(qs), dimnames = list(rownames=epsilons,colnames=qs))

## RUN
count <- 1
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      #for (l in 1:length(phis)) {
      for (m in 1:length(t1.changes)) {
        for (a in 1:length(gamma1s)) {
          for (b in 1:length(discount1s)) {
            params <- list( q=qs[i], 
                            epsilon=epsilons[j], 
                            db=dbs[k],
                            t1.response=ifelse(t1.changes[m]>=x$Tau,'NO',
                                               ifelse(gamma1s[a]>0,'CSR',
                                                      ifelse(discount1s[b]>0,'PRICE','?'))))
            print(sapply(names(params),function(name)sprintf("%s",params[name])))
            x$q <- qs[i]
            x$epsilon <- epsilons[j]
            x$db1 <- x$db2 <- dbs[k]
            x$gamma1 <- gamma1s[a]       ## RESPONSE ACTION: CSR
            x$discount1 <- discount1s[b] ## RESPONSE ACTION: PRICE
            x$t1.change <- t1.changes[m] ## RESPONSE TIMING
            index <- sprintf("q%s_epsilon%s_db%s_gamma1_%s_t1change%s_discount1_%s",
                             x$q, x$epsilon, x$db1, x$gamma1, x$t1.change, x$discount1)
            l.list[[index]] <- playCsrBayesGame(x, verbose = F)  
            l.list[[index]]$params <- params
            ##
            z[j,i] <- l.list[[index]]$B$B1[x$Tau] / sum(l.list[[index]]$B[x$Tau, ])
            count <- count+1            
          }
        }

      }
      #}
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file
# ## load saved binary file as follows:
# ## load('filename.RData')

idstr <- sprintf('_uber_l_list_facet_plot_CSR_NoResponse_T_%s_w_%s_J1_%s_J2_%s_gamma1_%s_discount1_%s_gn_%s_t1_%s',
                    x$Tau,x$omega,x$J1.0,x$J2.0,x$gamma1,x$discount1,n,paste(t1.changes,collapse="-")) 
saveRDS(l.list, file=sprintf("%s.rds",idstr))

## CREATE PLOT COORDS
tz <- t(z)
xx <- qs
yy <- epsilons

## SAVE PLOT
file.name <- sprintf('img/equilib_demand_share_%s.png', idstr)
png(file.name, height = 6, width = 7, units = 'in', res=250)
par(mar=c(4.2,4.2,3,1))
filled.contour(xx, yy, tz,  
               col=colorRamps::matlab.like2(21),
               xlab=expression('Buyer CSR Preference '~(q)),
               ylab=expression('Network Effect '~(epsilon)),
               main=sprintf("Platform 1 Equilibrium Demand Share\nCSR Cost %.2f, Response Pd %s, Hedonic Value %.2f", x$phi1, t1.changes[1],x$omega))
dev.off()



#=====================================================================
#----------------- 2. PRICE RESPONSE SIMULATION ----------------------
#---------------------------------------------------------------------

##----RESPONSE ACTION (platform 2) ------------
gamma1s    <- 0   # CSR action (0 < gamma1 < (1-c1)) #where c is marginal cost
discount1s <- .2  # PRICE  (0 < discount1 <= epsilon)
t1.changes <- t2r # x$Tau  #  c(30, 90, 270) # c(20, 80, 320) 
##-----FIXED PARAMS-----------------------------
dbs <-  0.1 # sellers churn # c(0.05, 0.1, 0.2) #c(.05,.5)
##----CONTINGENCIES TO TEST---------------------
qs <- round(seq(0,1, length.out = n), 2)  # seq(0,1,.1) ##c(0,.02,.05,.1,.2,.3,.4) # c(0,.2,.4,.6,.8,1)  #seq(0,1,.1)
epsilons <- round(seq(0.2,2, length.out = n), 2)  #
##----init data objects-------------------------
l.list <- list()  
z <- matrix(NA, nrow = length(epsilons), ncol=length(qs), dimnames = list(rownames=epsilons,colnames=qs))

## RUN
count <- 1
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      #for (l in 1:length(phis)) {
      for (m in 1:length(t1.changes)) {
        for (a in 1:length(gamma1s)) {
          for (b in 1:length(discount1s)) {
            params <- list( q=qs[i], 
                            epsilon=epsilons[j], 
                            db=dbs[k],
                            t1.response=ifelse(t1.changes[m]>=x$Tau,'NO',
                                               ifelse(gamma1s[a]>0,'CSR',
                                                      ifelse(discount1s[b]>0,'PRICE','?'))))
            print(sapply(names(params),function(name)sprintf("%s",params[name])))
            x$q <- qs[i]
            x$epsilon <- epsilons[j]
            x$db1 <- x$db2 <- dbs[k]
            x$gamma1 <- gamma1s[a]       ## RESPONSE ACTION: CSR
            x$discount1 <- discount1s[b] ## RESPONSE ACTION: PRICE
            x$t1.change <- t1.changes[m] ## RESPONSE TIMING
            index <- sprintf("q%s_epsilon%s_db%s_gamma1_%s_t1change%s_discount1_%s",
                             x$q, x$epsilon, x$db1, x$gamma1, x$t1.change, x$discount1)
            l.list[[index]] <- playCsrBayesGame(x, verbose = F)  
            l.list[[index]]$params <- params
            ##
            z[j,i] <- l.list[[index]]$B$B1[x$Tau] / sum(l.list[[index]]$B[x$Tau, ])
            count <- count+1            
          }
        }

      }
      #}
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file
# ## load saved binary file as follows:
# ## load('filename.RData')
idstr <- sprintf('_uber_l_list_facet_plot_CSR_PRICE_T_%s_w_%s_J1_%s_J2_%s_gamma1_%s_discount1_%s_gn_%s_t1_%s',
                 x$Tau,x$omega,x$J1.0,x$J2.0,x$gamma1,x$discount1,n,paste(t1.changes,collapse="-")) 
saveRDS(l.list, file=sprintf("%s.rds",idstr))

## CREATE PLOT COORDS
tz <- t(z)
xx <- qs
yy <- epsilons

## SAVE PLOT
file.name <- sprintf('img/equilib_demand_share_%s.png', idstr)
png(file.name, height = 6, width = 7, units = 'in', res=250)
par(mar=c(4.2,4.2,3,1))
filled.contour(xx, yy, tz,  
               col=colorRamps::matlab.like2(21),
               xlab=expression('Buyer CSR Preference '~(q)),
               ylab=expression('Network Effect '~(epsilon)),
               main=sprintf("Platform 1 Equilibrium Demand Share\nCSR Cost %.2f, Response Pd %s, Hedonic Value %.2f", x$phi1, t1.changes[1],x$omega))
dev.off()




#=====================================================================
#----------------- 3. CSR RESPONSE SIMULATION ----------------------
#---------------------------------------------------------------------

##----RESPONSE ACTION (platform 2) ------------
gamma1s    <- .2  # CSR action (0 < gamma1 < (1-c1)) #where c is marginal cost
discount1s <- 0   # PRICE  (0 < discount1 <= epsilon)
t1.changes <- t2r # x$Tau  #  c(30, 90, 270) # c(20, 80, 320) 
##-----FIXED PARAMS-----------------------------
dbs <-  0.1 # sellers churn # c(0.05, 0.1, 0.2) #c(.05,.5)
##----CONTINGENCIES TO TEST---------------------
qs <- round(seq(0,1, length.out = n), 2)  # seq(0,1,.1) ##c(0,.02,.05,.1,.2,.3,.4) # c(0,.2,.4,.6,.8,1)  #seq(0,1,.1)
epsilons <- round(seq(0.2,2, length.out = n), 2)  #
##----init data objects-------------------------
l.list <- list()  
z <- matrix(NA, nrow = length(epsilons), ncol=length(qs), dimnames = list(rownames=epsilons,colnames=qs))

## RUN
count <- 1
for (i in 1:length(qs)) {
  for (j in 1:length(epsilons)) {
    for (k in 1:length(dbs)) {
      #for (l in 1:length(phis)) {
      for (m in 1:length(t1.changes)) {
        for (a in 1:length(gamma1s)) {
          for (b in 1:length(discount1s)) {
            params <- list( q=qs[i], 
                            epsilon=epsilons[j], 
                            db=dbs[k],
                            t1.response=ifelse(t1.changes[m]>=x$Tau,'NO',
                                               ifelse(gamma1s[a]>0,'CSR',
                                                      ifelse(discount1s[b]>0,'PRICE','?'))))
            print(sapply(names(params),function(name)sprintf("%s",params[name])))
            x$q <- qs[i]
            x$epsilon <- epsilons[j]
            x$db1 <- x$db2 <- dbs[k]
            x$gamma1 <- gamma1s[a]       ## RESPONSE ACTION: CSR
            x$discount1 <- discount1s[b] ## RESPONSE ACTION: PRICE
            x$t1.change <- t1.changes[m] ## RESPONSE TIMING
            index <- sprintf("q%s_epsilon%s_db%s_gamma1_%s_t1change%s_discount1_%s",
                             x$q, x$epsilon, x$db1, x$gamma1, x$t1.change, x$discount1)
            l.list[[index]] <- playCsrBayesGame(x, verbose = F)  
            l.list[[index]]$params <- params
            ##
            z[j,i] <- l.list[[index]]$B$B1[x$Tau] / sum(l.list[[index]]$B[x$Tau, ])
            count <- count+1            
          }
        }

      }
      #}
    }
  }
}

# ## UNCOMMENT to SAVE binary (.RData) data file
# ## load saved binary file as follows:
# ## load('filename.RData')
idstr <- sprintf('_uber_l_list_facet_plot_CSR_CSR_T_%s_w_%s_J1_%s_J2_%s_gamma1_%s_discount1_%s_gn_%s_t1_%s',
                 x$Tau,x$omega,x$J1.0,x$J2.0,x$gamma1,x$discount1,n,paste(t1.changes,collapse="-")) 
saveRDS(l.list, file=sprintf("%s.rds",idstr))

## CREATE PLOT COORDS
tz <- t(z)
xx <- qs
yy <- epsilons

## SAVE PLOT
file.name <- sprintf('img/equilib_demand_share_%s.png', idstr)
png(file.name, height = 6, width = 7, units = 'in', res=250)
par(mar=c(4.2,4.2,3,1))
filled.contour(xx, yy, tz,  
               col=colorRamps::matlab.like2(21),
               xlab=expression('Buyer CSR Preference '~(q)),
               ylab=expression('Network Effect '~(epsilon)),
               main=sprintf("Platform 1 Equilibrium Demand Share\nCSR Cost %.2f, Response Pd %s, Hedonic Value %.2f", x$phi1, t1.changes[1],x$omega))
dev.off()











##_--------------------------------------------------------------
# kde <- kde2d(xx, yy, n=n)
# k <- n
# my.cols <- rev(brewer.pal(k, "RdYlBu"))
# padx <- .1; pady <- .1
# xlim <- c(min(xx)/(1+padx), max(xx)*(1+padx))
# ylim <- c(min(yy)/(1+pady), max(yy)*(1+pady))
# 
# contour(z, xlim = xlim, ylim=ylim)
# 
# image(z, ylim=ylim, xlim=xlim)
# 
# contour(kde, drawlabels=T, nlevels=k, col=my.cols, add=F,lwd=2)
# 
# # contour(x,y,z, add=F)
# 
# image(z, xlim=xlim, ylim=ylim)

##--------- PLOT BUYER SHARE -------------------------

## choose params to display
db_i <- dbs   # "0.05", "0.1", "0.2"
phi_i <- phis # "0.4" #c("0.01","0.1","1")
q_i <- qs # c('0', '0.2', '0.4', '0.6', '0.8', '1')  # .4
t1.change_i <- t1.changes


## subset data by chosen params
df <- getBasePlotDf(l.list, id.vars=c(names(l.list[[1]]$params), 'period') )
df <- subset(df, subset=(db %in% db_i & phi %in% phi_i & t1.change %in% t1.change_i))
df <- subset(df, subset=(as.character(q) %in% q_i))  ## couldn't filter 0.6 as float ?? need to use characters

## set numerics to factors for plotting
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
xintercepts1 <- l.list[[1]]$t1.change+1
colourCount <- length(unique(df$q))
db_num <- as.numeric(db_i)
db_lab <- ifelse(db_num < 0.1, 'Low',ifelse(db_num < 0.2, 'Moderate','High'))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# make labeller
labs <- c('Low','Moderate','High')
labels.epsilon = sapply(seq_along(epsilons), function(x)sprintf("%s Network Effect = %3s",labs[x],epsilons[x]))
labels.phis = sapply(seq_along(phi_i), function(x)sprintf("%s CSR Premium = %3s",labs[x],phi_i[x]))
names(labels.epsilon) <- epsilons
names(labels.phis) <- phi_i
labels <- c(labels.epsilon, labels.phis)

# legend guide
legend.guide <- guide_legend(title="Hedonic\nProportion\n(q)")

##----------- No title figure for article ----------------------------
## GGPLOT object
gg1 <- ggplot(aes(x=period, y=value, colour=q), data=df) + 
  geom_line(aes(colour=q,group=q,lty=q), lwd=1.1) + 
  facet_grid(phi ~ epsilon, labeller = as_labeller(labels)) + 
  geom_point(aes(pch=q), data=df.point)   +
  scale_color_manual(values=colorRamps::matlab.like(colourCount))  +
  geom_vline(xintercept=xintercepts2, lty=1) +
  geom_vline(xintercept=xintercepts1, lty=2) +
  scale_x_log10() + ylim(0,1) +
  ylab("Buyer Share") + xlab(expression('Time Period ('*log[10]~scale*')')) +
  guides(color=legend.guide, group=legend.guide,lty=legend.guide,pch=legend.guide) +
  theme_bw() 
gg1  ## display plot

## save plot
file.name <- sprintf("_uber_buyer_base_share_PRICE_6p5-10_not_log_J1_%s_J2_%s_t1_%s_t2_%s_T%s_db_%s_phi_%s_omega_%s.png", x$J1.0, x$J2.0,x$t1.change,x$t2.change,x$Tau,db_i,paste(phi_i,collapse="-"),x$omega)
ggsave(file.name, gg1, height=6.5, width=10, units='in')
##------------------------------------------------------------

