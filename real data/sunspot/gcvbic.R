require ( plot3D )
require (fda )
require ( latex2exp )
require ( pracma )
require ( mgcv )

setwd("E:/paper/functional period/code")

source ('my.bibasis1.R')

TTrend <- function (s,t){
  #z <- 25*t*sin (2*pi*s)
  z <- 2*(s+4*t)^2-2*s^2-32/3-8*s
  #z<-20*sin(2*pi*t+s)
  return (z) }
TPeriod <- function (s,t,theta0){
  #z <- 3*s^2**cos(2*pi*t/theta0+1.5*pi)
  z<-20*cos(2*pi*t/theta0+10*s)
  #z<-20*cos(2*pi*t/theta0+2*s)
  #z<-10*cos(2*pi*t/theta0)*s^2
  return (z) }

# TPeriod1 <- function (s,t){
#   #z <- 3*s^2**cos(2*pi*t/theta0+1.5*pi)
#   z<-2*cos(2*pi*t/12+s)
#   return (z) }
Tkernel <- function (t,u, norm.op =0.5 ){
  C <- norm.op/ 0.746824
  z <- C * exp( -(t^2 + u^2) / 2 )
  return (z)}


#The functional white noise Wn(s) is simulated as a Brownian motion in [0, 1],
data.ind <- function (N=100 ,m=50 , TT =1,sigma=1) {
  #N = sample size
  #m = num of points in [0 ,1]
  data <- NULL
  tt= seq (1,m, by=TT)/m
  data <- brow (TT,m,N,sigma=sigma)
  return ( data ) }

brow <- function (TT,m,N, sigma =1) {
  #TT= [0, TT]
  h=1/m
  Aux = matrix ( rnorm (n=N*(m-1), sd= sqrt (h* sigma ) ), m-1,N )
  Bt= rbind ( rep (0,N),apply (Aux, 2, cumsum ) )
  return (Bt) }


functTSTrend <- function (N=100, m=50, norm.op =0.5,theta0=12,sigma=1, burn =50) {
  # N= sample size
  # m= num of points in [0 ,1]
  h=1/m
  s=seq (1,m, by=1)/m
  KK <- outer (s,s, function (t,s) Tkernel (t,s, norm.op = norm.op))
  fw <- data.ind (N=N+ burn +1,m=m,TT =1,sigma=sigma)
  
  X0=fw [ ,1]
  fw=fw [,-1]
  DatosFAR <- matrix (0, length (s),N+ burn )
  DatosFAR [ ,1] <- X0
  for (i in 2:(N+ burn )){
    DatosFAR [,i] <- (t(KK) %*% DatosFAR [,(i -1) ])*h + fw[,i]
  }
  DatosFAR <- DatosFAR [,( burn +1) :(N+ burn )] # functional TS
  # simulating trend
  t_m=seq (1,N, by=1)
  t <- t_m/N
  trendM <- outer (s,t, function (s,t) TTrend (s,t))
  periodM <- outer (s,t_m, function (s,t_m) TPeriod(s,t_m,theta0))
  DatosM <- DatosFAR +trendM+periodM
  return ( list ( DataM =DatosM, TrendM =trendM, PeriodM=periodM,DataX = DatosFAR, KernelM =KK))
}
calx<-function(thet,T){
  x<-diag(thet)
  q<-ceiling(T/thet)
  m<-x
  for(i in 1:q){
    X1<-rbind(m,x)
    m=X1
  }
  X<-X1[1:T,]
  if (thet==1) {X=matrix(X,T,1)}
  return(X)
}


simulation_gcv<-function(N,m,N_sim, norm.op,theta0,K11,K22,K33,c_a,sigma){
  gcv_bic=matrix(0,(K11-3)*(K22-3)*(K33-4),2)
  kseq=matrix(0,(K11-3)*(K22-3)*(K33-4),3)
  s <- seq (1,m, by=1)/m
  t <- seq (1,N, by=1)/N
  t_m <- seq (1,N, by=1)
  # create basis for each coordinates
  
  Data <- functTSTrend (N,m, norm.op,theta0=theta0,sigma=sigma)
  YY <- Data$DataM
  #
  ymat<-matrix(YY, m*N, 1)
  
  
  X=calx(theta0,N)
  j_k=1
  for (K1 in 4:K11){
    for (K2 in 4:K22){
      for (K3 in 5:K33){
        kseq[j_k,1]=K1
        kseq[j_k,2]=K2
        kseq[j_k,3]=K3
  bases <- create.bspline.basis ( rangeval = c(0,1), nbasis = K1, norder = 4)
  basem <- create.bspline.basis ( rangeval = c(0,1), nbasis = K3, norder = 4)
  baset <- create.bspline.basis ( rangeval = c(0,1), nbasis = K2, norder = 4 )
  fdPs <- fdPar (bases, 2, lambda = 0 )
  fdPm <- fdPar (basem, 2, lambda = 0 )
  fdPt <- fdPar (baset, 2, lambda = 0 )
  # fdobjs  = fdPs$fd
  # sbasis  = fdobjs$basis
  # fdobjt  = fdPt$fd
  # tbasis  = fdobjt$basis
  sbasismat = eval.basis(s, bases)
  mbasismat = eval.basis(s, basem)
  tbasismat = eval.basis(t, baset)
  basismat  = kronecker(tbasismat,sbasismat)
  

    est.hat <- wsx.bibasis2 (K1,K2,K3,fdPs,fdPt, m,N, ymat,mbasismat, basismat,X)
    
    #Trend.hat.Mat <- est.hat$ghat_mat
    Period.hat.Mat <- est.hat$mhat_mat
    estgcv=est.hat$gcv
    estsse=est.hat$SSE
    bic_c= log(estsse)+log(m*N)/(2*m*N)*((K1)*(K2)+(K3*theta0))
    gcv_bic[j_k,1]=estgcv
    gcv_bic[j_k,2]=bic_c
    j_k=j_k+1
    
    }
  }
}
  #rate=o_sim/N_sim
  result0=list(gcv_bic=gcv_bic,kseq=kseq)
  return(result0)
  }
# plot
set.seed (1985)
N=70
m=50
N_sim=500
norm.op=0.5
sigma=1
theta0=12
K11=8
K22=8
K33=13
c_a=1
s <- seq (1,m, by=1)/m
t_m=seq (1,N, by=1)
t <- t_m/N
trendM <- outer (s,t, function (s,t) TTrend (s,t))
periodM <- outer (s,t_m, function (s,t_m) TPeriod(s,t_m,theta0))
s_time=proc.time()
resultgcv<-simulation_gcv(N,m,N_sim, norm.op,theta0,K11,K22,K33,c_a,sigma)
e_time=proc.time()

c(e_time-s_time)

c(which.min(resultgcv$gcv_bic[,1]),which.min(resultgcv$gcv_bic[,2]))

resultgcv$kseq[which.min(resultgcv$gcv_bic[,1]),]

resultgcv$kseq[which.min(resultgcv$gcv_bic[,2]),]

bicc=as.matrix(resultgcv$gcv_bic[,2],1)
gcvv=as.matrix(resultgcv$gcv_bic[,1],1)
ind_k=seq(1,length(bicc),by=1)
windows(width=14,height=10)
plot(ind_k,bicc,'l')

windows(width=14,height=10)
plot(ind_k,gcvv,'l')

resultgcv$kseq[which.min(resultgcv$gcv_bic[1:20,1]),]

resultgcv$kseq[which.min(resultgcv$gcv_bic[1:20,2]),]



