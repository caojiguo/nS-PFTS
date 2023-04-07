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


functTSTrend <- function (N=100, m=50, norm.op =0.5,theta0=12,sigma=1,burn =50) {
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


simulation<-function(N,m,N_sim, norm.op,theta0,K1,K2,K3,c_a,sigma){
etheta<-matrix(0,N_sim,3)
est_trend<-matrix(0,N_sim,N*m)
est_period<-matrix(0,N_sim,N*m)
ymat_sim<-matrix(0,N_sim,N*m)
o_sim=0
s <- seq (1,m, by=1)/m
t <- seq (1,N, by=1)/N
t_m <- seq (1,N, by=1)
# auxs <- rowMeans (YY)
# lambdas <- gam ( auxs ~te(s, bs='cr'),method ='REML')$sp
# auxt <- colMeans (YY)
# lambdat <- gam ( auxt ~te(t, bs='cr'),method ='REML')$sp
# sp.f=c( lambdas , lambdat )

# create basis for each coordinates

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
sk<-tk<-seq (1,200, by=1)/200
trend_eval=outer (sk,tk, function (sk,tk) TTrend (sk,tk))
for (j in 1:N_sim){
Data <- functTSTrend (N,m, norm.op,theta0,sigma)
YY <- Data$DataM
#
ymat<-matrix(YY, m*N, 1)
ymat_sim[j,]=ymat

#lse
Thm_0=min(40,max(24,N/2))
pseq<-seq(1,Thm_0,1)
Q2<-rep(0,Thm_0)
k<-1
for(thetai in pseq){
  try({
Xi=calx(thetai,N)
mySSE  = wsx.bibasis(ymat,mbasismat, basismat,Xi)})
Q2[k]<-mySSE
k<-k+1
  }
theta0ls = pseq[which.min(Q2)]
s1=min(Q2)/(m*N)

#bspline penalty

Q<-rep(0,Thm_0)
pseq<-seq(1,Thm_0,1)
k<-1
for(thetai in pseq){
  Q[k]=Q2[k]+c_a*(m*N)^(1/6)*log(N)*(thetai)*s1*m

  k=k+1
}
theta1 = pseq[which.min(Q)]
# pdf("QCcriterion5.pdf")
# #windows(width=14,height=10)
# plot(pseq,Q,"l",ylab='Q',xlab=expression(paste(theta)),cex.axis=1.5,cex.lab=1.5,lwd=2)
# dev.off()

# Trend, period estimation with estimated period
X=calx(theta1,N)
est.hat <- wsx.bibasis2 (K1,K2,K3,fdPs,fdPt, m,N, ymat,mbasismat, basismat,X)
  
#Trend.hat.Mat <- est.hat$ghat_mat
Period.hat.Mat <- est.hat$mhat_mat
#estgcv=est.hat$gcv

Trend.hat.Mat2 <- eval.bifd (sk,tk, est.hat$bifdobj )
rase_trend=as.vector(trend_eval-Trend.hat.Mat2)
rase2_trend=sqrt(mean(rase_trend^2))

rase_m=as.vector(Data$PeriodM-Period.hat.Mat)
rase2_m=sqrt(mean(rase_m^2))
etheta[j,1]=theta1
etheta[j,2]=rase2_trend
etheta[j,3]=rase2_m
# etheta[j,2]=ise2_trend
# etheta[j,3]=ise2_m
est_trend[j,]=est.hat$ghat
est_period[j,]=est.hat$mhat
if (theta1==theta0) o_sim<-o_sim+1

}
rate=o_sim/N_sim
result0=list(rate=rate,etheta=etheta,est_trend=est_trend,est_period=est_period,ymat_sim=ymat_sim)
return(result0)}
# plot
set.seed (1985)
N=70
m=24
N_sim=500
norm.op=0.3
sigma=1
theta0=12
K1=4
K2=5
K3=10
c_a=1
s <- seq (1,m, by=1)/m
t_m=seq (1,N, by=1)
t <- t_m/N
trendM <- outer (s,t, function (s,t) TTrend (s,t))
periodM <- outer (s,t_m, function (s,t_m) TPeriod(s,t_m,theta0))
s_time=proc.time()
result0<-simulation(N,m,N_sim, norm.op,theta0,K1,K2,K3,c_a,sigma)
e_time=proc.time()

c(e_time-s_time,result0$rate)

trendhat=colMeans(result0$est_trend)
periodhat=colMeans(result0$est_period)

Period.hat.Mat =matrix(periodhat,m,N)
Trend.hat.Mat=matrix(trendhat,m,N)

perc=function(x){quantile(x,c(0.025,0.975))}
mg2=apply(result0$est_trend,2,perc)
mp2=apply(result0$est_period,2,perc)
Trend.hat.Matl=matrix(mg2[1,],m,N)
Trend.hat.Matu=matrix(mg2[2,],m,N)
Period.hat.Matl=matrix(mp2[1,],m,N)
Period.hat.Matu=matrix(mp2[2,],m,N)

mrase2tr=mean(result0$etheta[,2])
mrase2pe=mean(result0$etheta[,3])
sdrase2tr=sd(result0$etheta[,2])
sdrase2pe=sd(result0$etheta[,3])

c(mrase2tr,sdrase2tr,mrase2pe,sdrase2pe)

Mm <- mesh (s,t)
mx <- Mm$x; my <- Mm$y
windows(width=14,height=10)
zg1=min(min(trendM),min(Trend.hat.Mat))
zg2=max(max(trendM),max(Trend.hat.Mat))
surf3D ( mx,my, Trend.hat.Mat, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7 
         #,plot=FALSE
         )

surf3D ( mx,my, Trend.hat.Matl, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7,
         add=TRUE,plot=FALSE)
#windows(width=14,height=10)
surf3D ( mx,my, Trend.hat.Matu, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7,
         add=TRUE,plot=TRUE)

#windows(width=14,height=10)
surf3D ( mx,my, trendM, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),   theta = 45, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('True$T(s,t/T)'),scale =TRUE , expand =.7 )

Mm1 <- mesh (s,t_m[1:theta0])
mx <- Mm1$x; my <- Mm1$y
zm1=min(min(periodM[,1:theta0]),min(Period.hat.Mat[,1:theta0]))
zm2=max(max(periodM[,1:theta0]),max(Period.hat.Mat[,1:theta0]))
windows(width=14,height=10)
surf3D ( mx,my, periodM[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
      zlim=c(zm1,zm2),   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('True$m(s,t)'),scale =TRUE , expand =.7 )

windows(width=14,height=10)
surf3D ( mx,my, Period.hat.Mat[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('\\hat{m}(s,t)'),scale =TRUE , expand =.7 
         #,plot=FALSE
         )

Mm <- mesh (s,t)
mx <- Mm$x; my <- Mm$y
windows(width=14,height=10)
zg1=min(min(trendM),min(Trend.hat.Mat))
zg2=max(max(trendM),max(Trend.hat.Mat))
surf3D ( mx,my, Trend.hat.Mat, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7 
         ,plot=FALSE
)

surf3D ( mx,my, Trend.hat.Matl, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7,
         add=TRUE,plot=FALSE)
#windows(width=14,height=10)
surf3D ( mx,my, Trend.hat.Matu, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),      theta = 45, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         main = TeX ('$\\hat{g}(s,t/T)'),scale =TRUE, expand =.7,
         add=TRUE,plot=TRUE)

Mm1 <- mesh (s,t_m[1:theta0])
mx <- Mm1$x; my <- Mm1$y
zm1=min(min(periodM[,1:theta0]),min(Period.hat.Mat[,1:theta0]))
zm2=max(max(periodM[,1:theta0]),max(Period.hat.Mat[,1:theta0]))
windows(width=14,height=10)
surf3D ( mx,my, periodM[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('True$m(s,t)'),scale =TRUE , expand =.7 )

windows(width=14,height=10)
surf3D ( mx,my, Period.hat.Mat[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('\\hat{m}(s,t)'),scale =TRUE , expand =.7 
         ,plot=FALSE)
surf3D ( mx,my, Period.hat.Matl[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('\\hat{m}(s,t)'),scale =TRUE , expand =.7 
         ,add=TRUE,plot=FALSE)
surf3D ( mx,my, Period.hat.Matu[,1:theta0], colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('\\hat{m}(s,t)'),scale =TRUE , expand =.7 
         ,add=TRUE,plot=TRUE)


Mmm <- mesh (s,t_m)
mx <- Mmm$x; my <- Mmm$y
zm1=min(min(periodM[,1:theta0]),min(Period.hat.Mat[,1:theta0]))
zm2=max(max(periodM[,1:theta0]),max(Period.hat.Mat[,1:theta0]))
windows(width=14,height=10)
surf3D ( mx,my, periodM, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  theta = 45, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('Truth$m(s,t)'),scale =TRUE , expand =.7 )
windows(width=14,height=10)
surf3D ( mx,my, Period.hat.Mat, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),   theta = 45, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         main = TeX ('\\hat{m}(s,t)'),scale =TRUE , expand =.7 )


##introduction toy example

TPeriod <- function (s,t,theta0){
  #z <- 3*s^2**cos(2*pi*t/theta0+1.5*pi)
  z<-sin(2*pi*t/theta0)*cos(2*pi*s)*t
  #z<-20*cos(2*pi*t/theta0+2*s)
  #z<-10*cos(2*pi*t/theta0)*s^2
  return (z) }

N=15
m=24
N_sim=500
norm.op=0.3
sigma=1
theta0=5
K1=4
K2=5
K3=10
c_a=1
s <- seq (1,m, by=1)/m
t_m=seq (1,N, by=1)
t_1=seq (1,theta0, by=1)
t <- t_m/N
#trendM <- outer (s,t, function (s,t) TTrend (s,t))
periodM1 <- outer (s,t_1, function (s,t_1) TPeriod(s,t_1,theta0))
periodM1=as.matrix(periodM1)
periodM=cbind(periodM1,periodM1,periodM1)
Mm1 <- mesh (s,t_m)
mx <- Mm1$x; my <- Mm1$y
zm1=min(periodM[,])
zm2=max(periodM[,])

pdf("toyfts.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, periodM, colkey = FALSE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), zlab='m(s,t)',  theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",scale =TRUE , 
         expand =.7,cex.axis=1.5,cex.lab=1.5)
# axis(1,seq(1,N,by=theta0),labels=seq(1,N,by=theta0))
dev.off()


pdf("toyfts2.pdf", width=14, height=7) 
#windows(width=14,height=7)
surf3D ( mx,my, periodM, colkey = FALSE, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2), zlab='m(s,t)',  theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",scale =TRUE , 
         expand =.7,cex.axis=1.5,cex.lab=1.5)
# axis(1,seq(1,N,by=theta0),labels=seq(1,N,by=theta0))
dev.off()

periodm=colMeans(periodM)

pdf("toymt.pdf")
#windows(width=7,height=7)
par(mgp=c(2,1,0))
plot(t_m,periodm,'l',ylim=c(zm1,zm2),xlab ="t", ylab =expression(paste(bar(m),'(t)')),cex.axis=1.5,cex.lab=1.5)
dev.off()

mvector=as.vector(periodM)

ts= seq (1,m*N, by=1)


x00=seq(12,360, 24)
x000=c(0,x00)
pdf("toymvector.pdf")
#windows(width=7,height=7)
plot(ts,mvector,'l',ylim=c(zm1,zm2), ylab ="m",lwd=2,xaxt="n",cex.axis=1.5,cex.lab=1.5,ann = F)
abline(v=0,lty=4)
abline(v=m,lty=4)
abline(v=2*m,lty=4)
abline(v=3*m,lty=4)
abline(v=4*m,lty=4)
abline(v=5*m,lty=5,col=2)
abline(v=6*m,lty=4)
abline(v=7*m,lty=4)
abline(v=8*m,lty=4)
abline(v=9*m,lty=4)
abline(v=10*m,lty=5,col=2)
abline(v=11*m,lty=4)
abline(v=12*m,lty=4)
abline(v=13*m,lty=4)
abline(v=14*m,lty=4)
abline(v=15*m,lty=5,col=2)
axis(1, tck=0,at=x00,labels=seq(1,15, 1),cex.axis=1.5,cex.lab=1.5,lwd=1)
title( ylab = 'm',cex.axis=1.5,cex.lab=1.5)
dev.off()
# 

