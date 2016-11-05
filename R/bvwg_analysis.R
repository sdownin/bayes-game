##
#
#  Bayesian Vote with the Wallet Game Analysis
#
##
source(file.path(getwd(),'R','bvwg_r6_game.R'))      # {.GLOBAL} GAME


getModelStr <- function() {
  "model{
    for (i in 1:n) {
      z[i] ~ dbern(q)
      th1[i] <- p2*(v1 + omega*sig1*z[i])
      th2[i] <- p1*(v2 + omega*sig2*z[i])
      s[i] <- ((th1[i]/th2[i])*pow(J1, epsilon)) / ( (th1[i]/th2[i])*pow(J1, epsilon) + pow(J2, epsilon) )
      G[i] ~ dbinom( s[i], L )
    }
    epsilon ~ dgamma(shapet,ratet)
    q ~ dbeta(h1t,h2t)
  }"
}

# getModelStr <- function() {
#   "model{
#     for (i in 1:n) {
#       z[i] ~ dbern(q)
#       th1[i] <- p2*(v1 + omega*sig1*z[i])
#       th2[i] <- p1*(v2 + omega*sig2*z[i])
#       X[i] <- dbern( ((th1[i]/th2[i])*pow(J1, epsilon)) / ((th1[i]/th2[i])*pow(J1, epsilon)+pow(J2, epsilon)) )
#     }
#     epsilon ~ dgamma(shapet,ratet)
#     q ~ dbeta(h1t,h2t)
#   }"
# }

dl2 <- data.list

jg <- JAGS$new(getModelStr(), c('q','epsilon'), list(name='Test1',pds=1))
jg$run(data.list, period=1 )
jg$bivarPlot('q','epsilon')


