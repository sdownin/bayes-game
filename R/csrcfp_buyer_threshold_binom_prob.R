setwd('C:\\Users\\sdowning\\Google Drive\\PhD\\Dissertation\\5. platform differentiation\\bayes-game')

library(ggplot2)


size <- 80
a1 = .5
a2 = .5
Ns <- c(1,10, 100,1000,10000)
mus <- c(0.05, 0.25, 0.5,0.75,.95)
qstars <- c(0.05, 0.25, 0.5,0.75,.95)
m <- length(Ns)*length(mus)*length(qstars)*size
df <- data.frame(prob=rep(NA,m),mu=rep(NA,m),
                 qstar=rep(NA,m),N=rep(NA,m),
                 x_N=rep(NA,m))
i <- 1
for (N in Ns) {
  for (mu in mus) {
    for (qstar in qstars) {
      n_is <- round(seq(1,N,length.out =size))
      for (n_i in n_is) {
        x <- round(qstar*(n_i+a1+a2)-a1+1)
        df$prob[i] <- dbinom(x, N, prob = mu)
        df$mu[i] <- sprintf('Buyer Prob (mu)=%s',mu)
        df$qstar[i] <- sprintf('q*=%s',qstar)
        df$N[i] <- sprintf('N=%s',N)
        df$x_N[i] <- x/N
        i <- i+1
      }
    }
  }
}
dfsub <- df[ df$x_N>=0 & df$x_N<=1 , ]
ggplot(aes(x=x_N,y=prob, colour=N), data=dfsub) +  facet_grid(mu ~ qstar) + 
  #geom_bar(stat='identity',position = 'dodge') +
  geom_line(aes(colour=N)) + geom_point(aes(colour=N), lwd=1.2) +
  xlab('x/N = (q*(N+a1+a2)-a1+1)/N') +
  ylab(expression(Binom(x*';'*mu[k]*','*N)%prop%kappa*'*')) +
  theme_bw()
##
file.name <- sprintf('Ns%s_PMF_buyer_votes_for_kstar.png',paste(Ns,collapse = "-"))
ggsave(file.name,  height=9,width=11,units = 'in',dpi=200)




mutrue <- 0.005

size <- 400
a1 = .5
a2 = .5
Ns <- c(2e7,4e7,8e7)
mus <- c(0.004, 0.0045, 0.005,0.0055,0.006)
qstars <- c(0.005,0.007, 0.009, 0.011,0.013)
m <- length(Ns)*length(mus)*length(qstars)*size
df <- data.frame(prob=rep(NA,m),mu=rep(NA,m),
                 qstar=rep(NA,m),N=rep(NA,m),
                 x_N=rep(NA,m))
i <- 1
for (N in Ns) {
  for (mu in mus) {
    for (qstar in qstars) {
      n_is <- seq(0,N,length.out = size)
      for (n_i in n_is) {
        x <- round(qstar*(n_i+a1+a2)-a1+1)
        df$prob[i] <- dbinom(x, N, prob = mu)
        df$mu[i] <- sprintf('Buyer Prob (mu)=%s',mu)
        df$qstar[i] <- sprintf('q*=%s',qstar)
        df$N[i] <- sprintf('N=%s',N)
        df$x_N[i] <- x/N
        i <- i+1
      }
    }
  }
}
dfsub <- df[ df$x_N>=0 & df$x_N<=1 , ]
ggplot(aes(x=x_N,y=prob, colour=N), data=dfsub) +  facet_grid(mu ~ qstar) + 
  #geom_bar(stat='identity',position = 'dodge') +
  geom_line(aes(colour=N)) + geom_point(aes(colour=N), lwd=1.2) +
  xlab('x/N = (q*(N+a1+a2)-a1+1)/N') +
  ylab(expression(Binom(x*';'*mu[k]*','*N)%prop%kappa*'*')) +
  xlim(0.003,0.007) +
  theme_bw()
##
file.name <- sprintf('Ns%s_PMF_buyer_votes_for_kstar.png',paste(Ns,collapse = "-"))
ggsave(file.name,  height=9,width=11,units = 'in',dpi=200)








df <- data.frame(prob=NA,mu=NA,qstar=NA,N=NA)

Ns.vec <- c(1,1000,1000000)
for (Ns in Ns.vec) {
  mus <- c(0.05, 0.1, 0.2)
  qstars <- c(.2,.4,.6)
  par(mfrow=c(length(mus),length(qstars)), oma=c(0,0,2,0))
  for (N in Ns) {
    for (i in 1:length(mus)) {
      for (j in 1:length(qstars)) {
        mu <- mus[i]
        qstar <- qstars[j]
        # qstar = .1
        # mu <- .05
        # N <- 1000
        Ns <- seq(0,N,length.out = 1000)
        a1 = a2 = .5
        x <- round(qstar*(Ns+a1+a2)-a1+1)
        barplot(height = dbinom(x, N, prob = mu), 
                names.arg = round(x/N,2),
                xlab='X/N', ylim=c(0,.25),
                main=sprintf('mu = %.2f, q* = %.2f', mu,qstar),
                col='steelblue', border='steelblue' )   
      }
    }
  }; title(sprintf('N = %s', N), outer=T)
}




library(lattice)
qstar = .1
N <- 1000
Ns <- seq(0,N,length.out = 80)
a1 = a2 = .5
x <- round(qstar*(Ns+a1+a2)-a1+1)
mu <- .1
df <- data.frame(x=factor(x), N=Ns, prob=dbinom(x, N, prob = mu))
dotplot(prob ~ x, data=df)
