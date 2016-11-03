setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\csr_bayes_game')
#source(file.path(getwd(),'R','csr_bayes_game_functions.R'))
library(fields)
library(pracma)

par(xpd=F)

#-------------------------------------------------------------------
#  PLATFORM CSR COSTS 
#-------------------------------------------------------------------
psi.star <- function(v1,v2,p1,p2,J1,J2,B1,B2,dB1,dB2,q,omega,rhotil)
{
  A1 <- (1+(1-dB1)*B1)*(v1+q*omega)*p2*J1^rhotil
  A2 <- (v2+q*omega)*v2*p1*J2^rhotil
  x1 <- (v1*v2*p2*J1^rhotil + A2) / ((v2+q*omega)*v1*p2*J1^rhotil + A2)
  x2 <- (A1 + (1-dB1)*B1*(v2+q*omega)*p1*J2^rhotil) / (A1 + (1-dB1)*B1*v2*p1*J2^rhotil)
  cat(sprintf('\ngam.star.1 < gam.star.2: %s\n',x1*x2<1))
  return (x1*x2)
}

## SIZE (SELLERS) DISADVANTAGE----------------------------------------
## INDIRECT NETWORK EFFECTS
## PRICE MARKUP
J1.vec <- seq(10,200,length.out = 20)
rhotil.vec <- 2^seq(-2,3,length.out = 6)
J2 <- 100
pst <- c()
for(rv in rhotil.vec){
  tmp <- psi.star(v1=1,v2=1,p1=20,p2=10,J1=J1.vec,J2=J2,B1=500,B2=500,dB1=.2,dB2=.2,q=.5,omega=1.5,rhotil=rv)
  pst <- cbind(pst, tmp)
}
file.title <- 'psi_star_sizeratio_price_markUP_rho_moderator_4.png'
png(file.path(getwd(),'img',file.title), height = 5, width = 6, units = 'in', res=250)
  par(mar=c(4.5,5,3,1))
  matplot(J1.vec/J2, pst, type='o',pch=15:20, 
          main=expression('Platform 1 Price Markup'~(p[1,t]/p[2,t])=='2/1'),
          xlab=expression('Size Ratio'~(J[1,t]/J[2,t])),ylab=expression('Pure Strategy Ratio'~(psi[1]^'*'/psi[1]^'**')))
  legend(x='topright',legend=rhotil.vec, title=expression(tilde(rho)),lty=1:6,pch=15:20, col=1:6)
  abline(h=1);abline(v=1)
dev.off()
## PRICE MARKDOWN
J1.vec <- seq(10,200,length.out = 20)
rhotil.vec <- 2^seq(-2,3,length.out = 6)
J2 <- 100
pst <- c()
for(rv in rhotil.vec){
  tmp <- psi.star(v1=1,v2=1,p1=10,p2=10,J1=J1.vec,J2=J2,B1=500,B2=500,dB1=.2,dB2=.2,q=.5,omega=1.5,rhotil=rv)
  pst <- cbind(pst, tmp)
}
file.title <- 'psi_star_sizeratio_prices_equal_rho_moderator_4.png'
png(file.path(getwd(),'img',file.title), height = 5, width = 6, units = 'in', res=250)
  par(mar=c(4.5,5,3,1))
  matplot(J1.vec/J2, pst, type='o',pch=15:20, 
          main=expression('Platform Prices Equal'~(p[1,t]~'='~p[2,t])),
          xlab=expression('Size Ratio'~(J[1,t]/J[2,t])),ylab=expression('Pure Strategy Ratio'~(psi[1]^'*'/psi[1]^'**')))
  legend(x='topright',legend=rhotil.vec, title=expression(tilde(rho)),lty=1:6,pch=15:20, col=1:6)
  abline(h=1);abline(v=1)
dev.off()
#
## PSI & GAMMA
file.title <- 'psi_gamma_star_sizeratio_rho_moderator_4.png'
png(file.path(getwd(),'img',file.title), height = 5, width = 6, units = 'in', res=250)
  par(mar=c(4.5,5,1,1))
  matplot(J1.vec/J2, pst, type='o',pch=15:20, xlab=expression('Size Ratio'~(J[1,t]/J[2,t])),ylab=expression('Pure Strategy Ratios'~(psi[1]^'*'/psi[1]^'**')~'and'~(gamma[1]^'*'/gamma[1]^'**')))
  legend(x='topright',legend=rhotil.vec, title=expression(tilde(rho)),lty=1:6,pch=15:20, col=1:6)
  abline(h=1);abline(v=1)
dev.off()
##---------------------------------
## PRICE w/ SIZE DISADVANTAGE 
p1.vec <- seq(10,200,length.out = 20)
rhotil.vec <- 2^seq(-2,3,length.out = 6)
p2 <- 100
pst <- c()
for(rv in rhotil.vec){
  tmp <- psi.star(v1=1,v2=1,p1=p1.vec,p2=p2,J1=50,J2=100,B1=500,B2=500,dB1=.2,dB2=.2,q=.5,omega=1.5,rhotil=rv)
  pst <- cbind(pst, tmp)
}
file.title <- 'psi_star_priceratio_disadvan_rho_moderator_1.png'
png(file.path(getwd(),'img',file.title), height = 5, width = 6, units = 'in', res=250)
  par(mar=c(4.5,5,3,1))
  matplot(p1.vec/p2, pst, type='o',pch=15:20, 
          main=expression('Platform 1 Size Disadvantage'~(J[1,t]/J[2,t])=='1/2'),
          xlab=expression('Price Ratio'~(p[1,t]/p[2,t])),ylab=expression('Pure Strategy Ratios'~(psi[1]^'*'/psi[1]^'**')~'and'~(gamma[1]^'*'/gamma[1]^'**')))
  legend(x='bottomright',legend=rhotil.vec, title=expression(tilde(rho)),lty=1:6,pch=15:20, col=1:6)
  abline(h=1);abline(v=1)
dev.off()
## PRICE w/ SIZE ADVANTAGE 
p1.vec <- seq(10,200,length.out = 20)
rhotil.vec <- 2^seq(-2,3,length.out = 6)
p2 <- 100
pst <- c()
for(rv in rhotil.vec){
  tmp <- psi.star(v1=1,v2=1,p1=p1.vec,p2=p2,J1=100,J2=50,B1=500,B2=500,dB1=.2,dB2=.2,q=.5,omega=1.5,rhotil=rv)
  pst <- cbind(pst, tmp)
}
file.title <- 'psi_star_priceratio_advan_rho_moderator_1.png'
png(file.path(getwd(),'img',file.title), height = 5, width = 6, units = 'in', res=250)
par(mar=c(4.5,5,3,1))
matplot(p1.vec/p2, pst, type='o',pch=15:20, 
        main=expression('Platform 1 Size Advantage'~(J[1,t]/J[2,t])=='2/1'),
        xlab=expression('Price Ratio'~(p[1,t]/p[2,t])),ylab=expression('Pure Strategy Ratios'~(psi[1]^'*'/psi[1]^'**')~'and'~(gamma[1]^'*'/gamma[1]^'**')))
legend(x='bottomright',legend=rhotil.vec, title=expression(tilde(rho)),lty=1:6,pch=15:20, col=1:6)
abline(h=1);abline(v=1)
dev.off()
##------------------------------------

## RESPONSE STRENGTH TO CSR
J1.vec <- seq(10,200,length.out = 20)
omega.vec <- seq(.1,4,length.out = 6)
J2 <- 100
pst <- c()
for(om in omega.vec){
  tmp <- psi.star(v1=1,v2=1,p1=10,p2=10,J1=J1.vec,J2=J2,B1=500,B2=500,dB1=.2,dB2=.2,q=.5,omega=om,rhotil=1.05)
  pst <- cbind(pst, tmp)
}
matplot(J1.vec/J2, pst, type='o',pch=20, xlab=expression(omega),ylab=expression(psi[1]^'*'/psi[1]^'**'))
legend(x='topright',legend=omega.vec, title=expression(tilde(omega)),lty=1:6,pch=20, col=1:6)
abline(h=1);abline(v=1)
## BUYER BASE 
J1.vec <- seq(10,200,length.out = 20)
B1.vec <- seq(10,2000,length.out = 6)
B2 <- 500
J2 <- 100
pst <- c()
for(B1 in B1.vec){
  tmp <- psi.star(v1=1,v2=1,p1=10,p2=10,J1=J1.vec,J2=J2,B1=B1,B2=B2,dB1=.2,dB2=.2,q=.5,omega=1.5,rhotil=1.05)
  pst <- cbind(pst, tmp)
}
matplot(J1.vec/J2, pst, type='o',pch=20, xlab=expression(omega),ylab=expression(psi[1]^'*'/psi[1]^'**'))
legend(x='topright',legend=B1.vec/B2, title=expression(B[1]/B[2]),lty=1:6,pch=20, col=1:6)
abline(h=1);abline(v=1)


## CHURN DISADVANTAGE----------------------------------------------------
dB1.vec <- seq(.01,.8,length.out = 50)
dB2 <- .2
pst <- psi.star(v1=1,v2=1,p1=10,p2=10,J1=50,J2=60,B1=500,B2=500,dB1=dB1.vec,dB2=dB2,
                q=.5,omega=1,rhotil=1)
plot(dB1.vec/dB2, pst, type='o',pch=20, xlab=expression(delta[B[1]]/delta[B[2]]),ylab=expression(psi[1]^'*'/psi[1]^'**'), ylim=c(.95,1.05))
abline(h=1);abline(v=1)


## CHURN DISADVANTAGE-----------------------------------------------------
B1.vec <- seq(10,1000,length.out = 50)
B2 <- 500
pst <- psi.star(v1=1,v2=1,p1=10,p2=10,J1=50,J2=48,B1=B1.vec,B2=B2,dB1=.2,dB2=.2,
                q=.5,omega=1,rhotil=1)
plot(B1.vec/B2, pst, type='o',pch=20, xlab=expression(B[1]/B[2]),ylab=expression(psi[1]^'*'/psi[1]^'**'), ylim=c(.9,1.1))
abline(h=1);abline(v=1)







# BETA
beta <- function(v1,v2,p1,p2,omega,q){
  x <- omega*q
  out <- ((v1+x)*p2)/((v2+x)*p1)
  return(out)
}

n <- 300
v1 <- 1
v2 <- seq(10,.5,length.out = n)
omega <- seq(.1,30,length.out = n)
grid <- expand.grid(v2=v2,omega=omega)

beta.vec <- apply(X = grid, MARGIN = 1, FUN = function(row){
  out <- beta(v1=v1,v2=row['v2'],p1=10,p2=20,omega=row['omega'],q=1)
  #print(out)
  return(out)
})
beta.mat <- matrix(beta.vec, nrow = n, byrow=FALSE)

#file.title <- 'hedonic_value_function_baseline_quality_omega_2.png'
#(file.path(getwd(),'img',file.title),height=5.5,width=6,units = 'in',res=250)
image.plot(v1/v2,omega,beta.mat, 
           xlab=expression(Baseline~Quality~Ratio~(upsilon[1]/upsilon[2])), 
           ylab=expression(Expected~Hedonic~Quality~(omega*q)~Ln~Scale),
           main=expression(Platform~Total~Quality~Ratio~(theta[1]*p[2*t]/theta[2]*p[1*t])), 
           log='y', legend.line = 2)
#dev.off()

contour(beta.mat, add=TRUE); box()


sh <- function(J1,omega){
  x1 <- (1+omega)*1*J1^1.1
  x2 <- (1+omega)*10*100^1.1
  s <- x1 / (x1+x2)
  return(s)
}


s.beta <- function(J1,J2,v1=1,v2=1,p1=10,p2=10,omega=5,q=1){
  x <- omega*q
  b <- ((v1+x)*p2)/((v2+x)*p1)
  return(b*J1 / (b*J1 + J2))
}

plot(s.beta(100,10, omega=seq(.1,20,length.out = 100)))


file.title <- 'hedonic_value_function_sizes_J1J2_2.png'
png(file.path(getwd(),'img',file.title),height=5.5,width=6,units = 'in',res=250)
image.plot(v1/v2,omega,beta.mat, 
           xlab=expression(Baseline~Quality~Ratio~(upsilon[1]/upsilon[2])), 
           ylab=expression(Expected~Hedonic~Quality~(omega*q)~Ln~Scale),
           main=expression(Platform~Total~Quality~Ratio~(theta[1]*p[2*t]/theta[2]*p[1*t])), 
           log='y', legend.line = 2)
dev.off()



#------ VECTOR FIELD  ----------------------------------------

drb1 <- function(B1,B2, s1=.5, m=.05, d1=.1, d2=.1){
  s2 <- 1-s1
  RB1 <- B1/B2
  M <- d1*B1+d2*B2 + (B1+B2)*m
  return((M/B2)*(s1-s2*RB1)+(d2-d1)*(RB1))
}

n <- 1000
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1))
for (s1 in seq(.1,.7,length.out = 4)) {
  drb1.s <- function(B1,B2) drb1(B1,B2,s1=s1)
  ## VECTORS
  vectorfield(drb1.s,c(0,n),c(0,n),
              lwd=2,scale = 1,col='steelblue',
              xlab=expression(B[1]),ylab=expression(B[2])
  ); abline(h=0);abline(v=0)
  ## PATHS
  for (xs in seq(0,n,length.out = 10)) {
    sol <- rk4(drb1.s, 0, n, xs, n/2)
    lines(sol$x, sol$y, col="darkgreen")
  }
}
grid()



##--EXAMPLE-----------------------------------
f <- function(x, y) x^2 - y^2
xx <- c(-1, 1); yy <- c(-1, 1)

## Not run:
vectorfield(f, xx, yy, scale = 0.1)
for (xs in seq(-1, 1, by = 0.25)) {
  sol <- rk4(f, -1, 1, xs, 100)
  lines(sol$x, sol$y, col="darkgreen")
}
grid()## End(Not run)
