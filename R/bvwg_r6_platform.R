library(R6)

#----------------------
# x <- list(
#      J = c(10)
#   ,  B = c(100)
#   ,  v1= 1
#   , db1=.5  # 30% buy all (y/pk) goods from current platform 1; 70% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
#   , dj1=.1
#   , c1=.5       ## seller MARGINAL cost
#   , gamma1=.05  ## seller CSR cost
#   , w1=.02      ## Platform operator MARGINAL cost
#   , psi1=.02    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
#   , a1=1
#   , r1=.1
#   , sig1.fixed=NA
#   , t1.change=NA
# )
# ## MAIN GAME CALL
# ## SET STRATEGY
# x$t1.change <- 1
# x$t2.change <- 4
# x$sig1.fixed <- c(rep(0,x$t1.change),rep(1,x$Tau-x$t1.change))
# x$sig2.fixed <- c(rep(0,x$t2.change),rep(1,x$Tau-x$t2.change))
#----------------------

PLATFORM <- R6Class("PLATFORM", 
                    public = list(
                      name = NA,
                      l = list(
                        sig = c(0)   #0=No CSR
                        , J = c(10)       #sellers
                        , B = c(100)      #buyers
                        , v= 1            #utilitarian value
                        , db=.5  # 50% buy all (y/pk) goods from current platform 1; 50% defect to multihome buying s1*(y/p1) from Plat 1, s2*(y/p2) from Plat 2
                        , dj=.1
                        , c=.5       ## seller MARGINAL cost
                        , gamma=.05  ## seller CSR cost
                        , w=.02      ## Platform operator MARGINAL cost
                        , psi=.02    ## Platform operator CSR cost   moved --> function of (gamma, B, y, p1)
                        , a=1
                        , r=.1
                      ),
                      initialize = function(name=NA, l=list()) {
                        self$name <- name
                        self$add(l)
                      },
                      
                      add = function(l) {
                        self$l <- c(self$l, l)
                      },
                      
                      set = function(item, val) {
                        self$l[[item]] <- val
                      },
                      
                      trans = function() {
                        
                      },
                      
                      getPsi = function(gamma,y,p,B) {
                        return(gamma / ((y/p)*B))
                      },
                      
                      getB = function(s,m,b,d)  {
                        newB <- s*m + b*(1-d)
                        return(ifelse(newB > 0, newB, 0))
                      },
                      
                      getJ = function(y,rho,gamma,c,B,f,J,dj)  {
                        netMargProfit <- ((rho+1)/rho) - (gamma/(rho*c))
                        newJ <- y*netMargProfit*(B/f) + J*(1-dj)
                        return(ifelse(newJ > 0, newJ, 0))
                      },
                      
                      getG = function(s,L,M,seed=1111)  {
                        set.seed(seed)
                        size <- min(M, length(s))
                        s.sample <- sample(s,size,replace = F)
                        return( sapply(s.sample, function(s)rbinom(n=1, size = L, s)) )
    },
    
    getQty = function(y,p,netB,G) {
      return(netB*(y/p) + G)
    },
    
    getPi = function(r,omega,psi,rho,c,Q,O,y)  {
      netMargProfit <-  r - ( (omega+psi)/(rho*c) )
      return( netMargProfit*Q - (O/y) )
    },
  
    muBeta = function(a,b)  {
      return(a/(a+b))
    },
    
    varBeta = function(a,b)  {
      num <- a*b
      denom <- (a+b)^2 * (a+b+1)
      return(num/denom)
    },
    
    getPrice = function(k, x)  {
      rhoTilde <- 1 + x$rho
      if(k %in% c(1,'1'))
        return( rhoTilde*x$c1 )
      if(k %in% c(2,'2'))
        return( rhoTilde*x$c2 )
    },
    
    getMaxPsi = function(k, x) {
      rhoTilde <- 1 + x$rho
      if(k %in% c(1,'1'))
        return(x$r1*rhoTilde*x$c1 - x$w1)
      if(k %in% c(2,'2'))
        return(x$r2*rhoTilde*x$c2 - x$w2)
    }
  )
)

