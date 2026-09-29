## k inert factors, range = [0,1]
inert <- function(xx, p, k, seed = 0, sd) {
  val = 0
  set.seed(seed)
  for(i in 1:k) {
    a = rnorm(1, mean = 0, sd = sd)
    val = val + a * xx[p+i]
  }
  for(i in 1:k) {
    for(j in 1:p) {
      b = rnorm(1, mean = 0, sd = sd)
      val = val + b * xx[p+i] * xx[j]
    }
  }
  return(val)
}

# g-function
f1 <- function(xx, a = c(rep(0,2), rep(3,6))) {
  y <- 1
  for(i in 1:8) {
    y = y*(abs(4*xx[i] - 2) + a[i]) / (1+a[i])
  }
  return(y)
}
g_fun <- function(X) apply(X, 1, f1)


# levy function
f2 <- function(xx)
{
  xx <- scales::rescale(xx, to = c(-10,10), from = c(0,1))
  xx <- xx[1:8]
  d <- length(xx)
  w <- 1 + (xx - 1)/4
  
  term1 <- (sin(pi*w[1]))^2 
  term3 <- (w[d]-1)^2 * (1+1*(sin(2*pi*w[d]))^2)
  
  wi <- w[1:(d-1)]
  sum <- sum((wi-1)^2 * (1+10*(sin(pi*wi+1))^2))
  
  y <- term1 + sum + term3
  return(y)
}
levy <- function(X) apply(X,1,f2)


# ACKLEY FUNCTION
f3 <- function(xx, a=20, b=0.2, c=2*pi)
{
  xx <- xx[1:8]
  xx <- scales::rescale(xx, to = c(-32.768, 32.768), from = c(0,1))
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}
ackley <- function(X) apply(X,1,f3)


# borehole function
f4=function(x) {
  x <- x[1:8]
  lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
  upper<-c(0.15,50000,115600,1110,116,820,1680,12045)
  x=lower+x*(upper-lower)
  val=2*pi*x[3]*(x[4]-x[6])/(log(x[2]/x[1])*(1+2*x[7]*x[3]/(log(x[2]/x[1])*x[1]^2*x[8])+x[3]/x[5]))
  return(val)
}
borehole=function(X) apply(X,1,f4)


# Dette and Pepelyshev function
f5 <- function(x){
  for(k in 4:8){
    term4=sum(k*log(1+sum(x[3:k])))
  }
  val=4*(x[1]-2+8*x[2]-8*x[2]^2)^2+(3-4*x[2])^2+16*sqrt(x[3]+1)*(2*x[3]-1)^2+term4
  return(val)
}
Dette <- function(X) apply(X,1,f5)


# OTL CIRCUIT FUNCTION
f6 <- function(xx)
{
  xx <- xx[1:6]
  xx[1] <- scales::rescale(xx[1], to = c(50,150), from = c(0,1))
  xx[2] <- scales::rescale(xx[2], to = c(25,70), from = c(0,1))
  xx[3] <- scales::rescale(xx[3], to = c(0.5,3), from = c(0,1))
  xx[4] <- scales::rescale(xx[4], to = c(1.2,2.5), from = c(0,1))
  xx[5] <- scales::rescale(xx[5], to = c(0.25,1.2), from = c(0,1))
  xx[6] <- scales::rescale(xx[6], to = c(50,300), from = c(0,1))
  
  Rb1  <- xx[1]
  Rb2  <- xx[2]
  Rf   <- xx[3]
  Rc1  <- xx[4]
  Rc2  <- xx[5]
  beta <- xx[6]
  
  Vb1 <- 12*Rb2 / (Rb1+Rb2)
  term1a <- (Vb1+0.74) * beta * (Rc2+9)
  term1b <- beta*(Rc2+9) + Rf
  term1 <- term1a / term1b
  
  term2a <- 11.35 * Rf
  term2b <- beta*(Rc2+9) + Rf
  term2 <- term2a / term2b
  
  term3a <- 0.74 * Rf * beta * (Rc2+9)
  term3b <- (beta*(Rc2+9)+Rf) * Rc1
  term3 <- term3a / term3b
  
  Vm <- term1 + term2 + term3
  return(Vm)
}
otlcircuit <- function(X) apply(X,1,f6)


# PISTON FUNCTION
f7 <- function(xx)
{
  xx <- xx[1:7]
  xx[1] <- scales::rescale(xx[1], to = c(30, 60), from = c(0,1))
  xx[2] <- scales::rescale(xx[2], to = c(0.005, 0.02), from = c(0,1))
  xx[3] <- scales::rescale(xx[3], to = c(0.002, 0.01), from = c(0,1))
  xx[4] <- scales::rescale(xx[4], to = c(1000, 5000), from = c(0,1))
  xx[5] <- scales::rescale(xx[5], to = c(90000, 110000), from = c(0,1))
  xx[6] <- scales::rescale(xx[6], to = c(290, 296), from = c(0,1))
  xx[7] <- scales::rescale(xx[7], to = c(340, 360), from = c(0,1))
  M  <- xx[1]
  S  <- xx[2]
  V0 <- xx[3]
  k  <- xx[4]
  P0 <- xx[5]
  Ta <- xx[6]
  T0 <- xx[7]
  
  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3
  
  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)
  
  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))
  
  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}
piston <- function(X) apply(X,1,f7)


# ROBOT ARM FUNCTION
f8 <- function(xx)
{
  xx <- xx[1:8]
  xx[1] <- scales::rescale(xx[1], to = c(0, 2*pi), from = c(0, 1))
  xx[2] <- scales::rescale(xx[2], to = c(0, 2*pi), from = c(0, 1))
  xx[3] <- scales::rescale(xx[3], to = c(0, 2*pi), from = c(0, 1))
  xx[4] <- scales::rescale(xx[4], to = c(0, 2*pi), from = c(0, 1))
  
  theta <- xx[1:4]
  L     <- xx[5:8]
  
  thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
  thetamatlow <- thetamat
  thetamatlow[upper.tri(thetamatlow)] <- 0
  sumtheta <- rowSums(thetamatlow)
  
  u <- sum(L*cos(sumtheta))
  v <- sum(L*sin(sumtheta))
  
  y <- (u^2 + v^2)^(0.5)
  return(y)
}
robot <- function(X) apply(X,1,f8)