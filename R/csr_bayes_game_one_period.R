
source(file.path(getwd(),"R","csr_bayes_game_functions_2.R"))
library(reshape2)


qpost <- function(W,N=50,a1=0.5,a2=0.5,theta=0.5) ## W = sum of votes with wallet
{
  y2 <- W + a2
  y1 <- N - W + a1
  z <- seq(0,1,.001)
  qW <- dbeta(z,y2,y1)
  qA <- dbeta(z,y1,y2)
  # qWqA <- qW * theta + qA * (1-theta)
  #qAdqW <- qA - qW
  df <- data.frame(qW=qW
                    #, qA=qA
                    #, qWqA=qWqA
                    #, qAdqW=qAdqW
                    )
  return(df)
}

## ONLY ONE VOTE WITH THE WALLET
par(mfrow=c(1,2), mar=c(4.5,4.5,4.5,1.5), oma=c(0,0,2,0))
wvec <- c(0,1)
adjs <- c(.8,.2)
frames <- c('(a) No CSR','(b) CSR')
for (i in seq_along(wvec)) {
  # W <- 2
  W <- wvec[i]
  N <- 1
  y <- qpost(W, N, a1=.5,a2=.5)
  x <- seq(0,1,length.out = nrow(y))
  colors <- c(4,'darkgray')
  # mainstr <- sprintf('%s  N = %s, h = %s',frames[i],N,W)
  mainstr <- frames[i]
  matplot(x, y, main=mainstr,
          xlab = 'Proportion of Hedonic Buyers', 
          ylab = expression('Posterior Hedonic Prob. '*(q^post)),
          type='o', pch=c(NA), lwd=2:1, lty=1,col=colors)
  qstar <- 1/2
  abline(v=qstar, lty=2)
  mtext(text = expression(q^"*"), at=qstar, side = 3, cex=.9, padj = .0)
  qWx <- y$qW* x
  meanqW <- mean( qWx[qWx>-Inf&qWx<Inf], na.rm = T )
  cat(sprintf('\nN=%s, W=%s, mean(qW) = %.3f\n',N,W,meanqW))
  abline(v=meanqW, lty=4)
  mtext(text = expression(hat(q)==E[X]*(q^post)), 
        at=meanqW, side = 3, cex=.9, adj=adjs[i],  padj = .0)
  leg.loc <- ifelse(W > 0, 'topleft', 'topright')
  legend(leg.loc,legend=c(expression(q^post)
                          #,expression(q^A)
                          #, expression(E[sigma^b==(W*","*A)]*q),
                          #, expression(q^A -q^W)
                          #,sprintf('h = %s',W),sprintf('N = %s',N)
  ),
  title=sprintf('%s (%s/%s)',ifelse(W>0,'vote','abstain'),W,N),
  pch=c(NA),lwd=c(2),lty=c(1),col=c(4))
}; mtext("Best Response to Observed Buyers Voting with the Wallet", outer=T, cex=1.5)


## MANY VOTES WITH WALLET
png('best_response_to_observed_buying_with_wallet.png', height=4.5, width=8.8, units='in', res=300)
  par(mfrow=c(1,2), mar=c(4.5,4.5,4.5,1.5), oma=c(0,0,2,0))
  wvec <- c(1,7)
  adjs <- c(.4,.1)
  frames <- c('(a) No CSR','(b) CSR')
  for (i in seq_along(wvec)) {
    # W <- 2
    W <- wvec[i]
    N <- 10
    y <- qpost(W, N, a1=.5,a2=.5)
    x <- seq(0,1,length.out = nrow(y))
    colors <- c(4,'darkgray')
    # mainstr <- sprintf('%s  N = %s, h = %s',frames[i],N,W)
    mainstr <- frames[i]
    matplot(x, y, 
            main=mainstr,
            xlab = 'Proportion of Hedonic Buyers', 
            ylab = 'Density',
            type='o', pch=c(NA), lwd=2:1, lty=1,col=colors)
    qstar <- 1/2
    abline(v=qstar, lty=1)
    mtext(text = expression(q*"*"), at=qstar, side = 3, cex=.9, padj = .0)
    qWx <- y$qW* x
    meanqW <- mean( qWx[qWx>-Inf&qWx<Inf], na.rm = T )
    cat(sprintf('\nN=%s, W=%s, mean(qW) = %.3f\n',N,W,meanqW))
    abline(v=meanqW, lty=4)
    mtext(text = expression(hat(q)==E[X]*(q^post)), 
          at=meanqW, side = 3, cex=.9, adj=adjs[i],  padj = .0)
    polygon(x=x, y=y$qW, col=rgb(.2,.2,.8,.2))
    leg.loc <- ifelse(W/N < .5, 'topright', 'topleft')
    legend(leg.loc,legend=c(expression(q^post)
                            , expression(q*'*')
                            , expression(hat(q))
                            #,expression(q^A)
                            #, expression(E[sigma^b==(W*","*A)]*q),
                            #, expression(q^A -q^W)
                            #,sprintf('h = %s',W),sprintf('N = %s',N)
          ),
          title=sprintf('%s votes (%s/%s)',ifelse(W/N>.35,'Many','Few'),W,N),
          pch=c(NA),lwd=c(3,1,1),lty=c(1,1,4),col=c(4,1,1))
  }; #mtext("Best Response to Observed Buyers Voting with the Wallet", outer=T, cex=1.5)
dev.off()

##------------------------------------------------------
getPriceNox <- function(sig,epsilon,gamma,phi)
{
  return((1+epsilon)*(1 + sig*gamma + sig*phi))
}
getPriceFromMarkup <- function(sig,epsilon,gamma,phiMarkup, mc=1)
{
  phi <- phiMarkup*(mc + sig*gamma)
  return((1+epsilon)*(mc + sig*gamma + sig*phi))
}

zeta2 <- function(z,sig1, sig2=1,epsilon=1.1,omega=1,u1=1,u2=1,cost1=1,cost2=1,gamma1=.2,gamma2=.2,phi1=.1,phi2=.1)
{
  a2 <- (u2 + omega*z*sig2) * (cost1 + (gamma1 + phi1)*z*sig1 )
  a1 <- (u1 + omega*z*sig1) * (cost2 + (gamma2 + phi2)*z*sig2 )
  expo <- 1/epsilon
  return( (a2/a1)^expo )
}
J1fromJ2 <- function(J2,z,sig1, epsilon=1.1)
{
  return(J2 * zeta2(z, sig1, epsilon=epsilon))
}
##------------------------------------------------------
a1 <- function(sig1, sig2, Z, u1=1, u2=1, omega=1, J1=200, J2=50, epsilon=1.05, phiMarkup=.1)
{
  gamma1 <- gamma2 <- 0.2
  c1 <- c2 <- 1
  p2 <- getPriceFromMarkup(sig2, epsilon, gamma2, phiMarkup, c2)
  return((u1 + omega*Z*sig1)*p2*J1^epsilon)
}
a2 <- function(sig1, sig2, Z, u1=1, u2=1, omega=1, J1=200, J2=50, epsilon=1.05, phiMarkup=.1)
{
  gamma1 <- gamma2 <- 0.2
  c1 <- c2 <- 1
  p1 <- getPriceFromMarkup(sig1, epsilon, gamma1, phiMarkup, c1)
  return((u2 + omega*Z*sig2)*p1*J2^epsilon)
}


xi1 <- function(sig1)
{
  y <- (a1(sig1,1,0) + a2(sig1,1,0)) / (a1(sig1,1,1) + a2(sig1,1,1))
  z <- (a2(sig1,1,1) - a2(sig1,1,1)) / (a2(sig1,1,1)*a1(0,1,1) - a1(sig1,1,1)*a2(sig1,0,1))
  return(y * z)
}

xi1 <- function(sig1, ...)
{
  y <- (a1(sig1,1,0,...) + a2(sig1,1,0,...)) / (a1(sig1,1,1,...) + a2(sig1,1,1,...))
  z <- (a2(sig1,1,1,...) - a2(sig1,1,1,...)) / (a2(sig1,1,1,...)*a1(0,1,1,...) - a1(sig1,1,1,...)*a2(sig1,0,1,...))
  return(y * z)
}

phis <- seq(0,1,.01)

plot(xi1(1,phiMarkup=phis))



##-------------------- zeta --------------------------------
png("profit_thresholds-1.png", height=7, width=9, units='in', res=250)
  epsilons <- c(1/10, 1/2, 1,4)
  J2 <- 2^(0:14)
  par( mfrow=c(2,2), oma=c(0,0,1,0), mar=c(4.1,4.1,1,1))
  z <- 0:1
  sig1 <- 0:1
  parms <- list(c(0,0), c(1,0), c(1,1))
  df <- data.frame()
  for (i in 1:length(epsilons)) {
    epsilon <- epsilons[i]
    J1 <- sapply(parms, function(x){
      cat(x); cat('\n')
      z <- x[1];  sig1 <- x[2]
      return(J1fromJ2(J2,z,sig1, epsilon))
    }); colnames(J1) <- c('U', 'H_NO', 'H_CSR')
    ##
    lims <- c(1,10000)
    matplot(J2, J1, ylim=lims,xlim=lims, log='xy', type='o', pch=1:4, lty=1:4, col=c(1,2,4), 
            ylab=expression("Platform 1 Size ("~J[1]~")"), xlab=expression("Platform 2 Size ("~J[2]~")"))
    mtext(text=expression(epsilon=="      "))
    mtext(text=sprintf('     %s',round(epsilon,2)))
    legs <- sapply(parms,function(x){
      paste0(ifelse(x[1]==1,'H','U'),ifelse(x[1]==0,'',ifelse(x[2]==1,', CSR',', NO')))
    })
    legend('bottomright', legend=c(legs, expression(J[1]==J[2])), col=c(1,2,4,'gray'), lty=c(1:3,4), pch=c(1:3,NA))
    text(expression(Delta*pi["1"]==a["1"](sigma["1"])*y), x = 6, y=5000)
    text(expression(Delta*pi["1"]==0), x = 100, y=2)
    abline(a = 0, b=1, lty=4, col='gray')
    ##
    df.tmp <-  data.frame(J1)
    df.tmp$J2 <- J2
    df.tmp$epsilon <- epsilon
    df <- rbind(df, df.tmp)
  }; title(main=expression("Profit Thresholds ("~J[1]==zeta[2]*(sigma[1]*",Z")*J[2]~")"), outer = T)
dev.off()

df_long <- melt(df, id.vars = c('epsilon', 'J2'))

ggplot(aes(x=J2, y=value), data=df_long) + 
  scale_x_log10() + scale_y_log10() +
  geom_line(aes(colour=variable)) +   geom_point(aes(colour=variable)) + 
  facet_wrap(~epsilon) + theme_bw() + ylab(expression(J[1])) + xlab(expression(J[2]))

xyplot(value ~ J2 | factor(epsilon), groups=variable, data=df_long, 
       main="Profit Thresholds",auto.key = T, 
       type="o", pch=1:3, col=c(1,2,4),ylab=expression(J[1]), xlab=expression(J[2]),
       scales=list(x=list(log=2),y=list(log=2)))

##----------------- PRICE CONTOUR -------------------------
epsilon <- seq(0,2,.02)
markup <- seq(0,2,.02)
getPriceFromMarkup.mesh <- Vectorize(function(markup, ep){
  return(getPriceFromMarkup(sig = 1, epsilon = ep, gamma = 0.2, phiMarkup = markup))
})
p <- outer(epsilon, markup, getPriceFromMarkup.mesh)
rownames(p) <- epsilon
colnames(p) <- markup

image(epsilon, markup, p, ylab=expression("Markup ("~varphi~"/p("~varphi~"))"), xlab=expression("Indirect Network Effects ("~epsilon~")"))
k <- 7; my.cols <- rev(brewer.pal(k, "RdYlBu"))
contour(epsilon, markup, p,  drawlabels=T, nlevels=k, col=my.cols, add=T, lwd=2)

#-------------------- SHARE CONTOUR (for platform 2 (1-s)) ------------------------
##  s = f(q,epsilon) 
ep.vec <- seq(0,2,.02)
q.vec <- seq(0,1,.01)
t <- 1

df <- data.frame()
for (omega in c(1,10)) {
  # png(sprintf("demand_share_%s_q.png",ifelse(omega==1,"Low", "High")), height=6, width=9, units='in', res=250)
  par(mfrow=c(2,3), mar=c(4.1,4.1,1.5,1), oma=c(0,0,2,0))
  counter <- 
  for (j in c(5,125,250,500,10500)) {
    share.mesh <- Vectorize(function(q,eps) {
      gamma1 <- gamma2 <- 0.2
      phi1 <- phi2 <- 0
      sig1 <- 0
      sig2 <- 1
      p1 <- getPriceNox(sig1, eps, gamma1, phi1)  #(1+eps)*(1 + sig1*gamma1 + sig1*phi1)
      p2 <- getPriceNox(sig2, eps, gamma2, phi2)  #(1+eps)*(1 + sig2*gamma2 + sig2*phi2)
      return(share(p1=p1,p2=p2,v1=1,v2=1,sig1=sig1,sig2=sig2,J1=j,J2=250,omega=omega,z=q,epsilon=eps,k=1))
    })
    q <- q.vec; epsilon <- ep.vec
    s <- 1-outer(q.vec, ep.vec, share.mesh)
    rownames(s) <- q
    colnames(s) <- epsilon
    df_long <- melt(s, varnames = c('q','epsilon'), value.name = "share")
    df_long$omega <- omega
    df <- rbind(df, df_long)
    ##
    image(q,epsilon,s, ylab=expression("Indirect Network Effects ("~epsilon~")"), xlab=expression("Hedonic Demand Proportion ("~q~")"), main=sprintf("Seller Share = %.3f",250/(j+250)))
    k <- 7; my.cols <- rev(brewer.pal(k, "RdYlBu"))
    contour(q,epsilon,s , drawlabels=T, nlevels=k, col=my.cols, add=T, lwd=2)
    ##
    counter <- counter + 1
  }; title(main=sprintf("Demand Share with %s Hedonic Value",ifelse(omega>1,"High","Low")), outer=T)
  # dev.off()
}

# df <- subset(df, subset=(omega==1))
# g1 <- ggplot(df, aes(q, epsilon, z=share)) #+ facet_wrap(~ omega)
# g1 + geom_contour()

persp.withcol(q,epsilon,s,heat.colors,20, 
              theta=140, phi=25, expand=0.5,
              ltheta=120, ticktype='detailed', border=NA, lty=6)


##


#--------------- PROFIT CONTOUR -------------------------------------

epsilon <- c(1.1)
##
q <- seq(0,1,.01)
markup <- seq(0,2,.02)
getPriceFromMarkup.mesh <- Vectorize(function(markup, ep){
  return(getPriceFromMarkup(sig = 1, epsilon = ep, gamma = 0.2, phiMarkup = markup))
})
p <- outer(epsilon, markup, getPriceFromMarkup.mesh)
rownames(p) <- epsilon
colnames(p) <- markup





#-------------------------------------------------------------------





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
