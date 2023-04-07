require ( plot3D )
require (fda )
require ( latex2exp )
require ( pracma )
require ( mgcv )
require ( ggplot2)
setwd("E:/paper/functional period/code")
##use
source ('my.bibasis1.R')
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
setwd("E:/paper/functional period/data")

mdata <- read.csv("hkair.csv")
TA=as.matrix(mdata[1:761,13])
ymat<-as.matrix(mdata[,7])
#ymat<-as.matrix(mdata[,8])
m<-24
N<-length(ymat)/m
YY<-matrix(ymat,m,N)

s <- seq (1,m, by=1)/m
t <- seq (1,N, by=1)/N
t_m <- seq (1,N, by=1)
K1=6
K3=9
K2=10
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
#lse
Thm_0=min(round(N/2),100)
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
c_a=1
Q<-rep(0,Thm_0)
pseq<-seq(1,Thm_0,1)
k<-1
for(thetai in pseq){
  Q[k]=Q2[k]+c_a*log(N*m)*(thetai)*s1*(N*m)^0.4/11
  k=k+1
}
theta1 = pseq[which.min(Q)]
pdf("Qhknox.pdf")
#windows(width=7,height=7)
plot(pseq,Q,"l",ylab='Q',xlab=expression(paste(theta)),cex.axis=1.5,cex.lab=1.5,lwd=2)
axis(1,theta1,labels=theta1,cex.axis=1.5,cex.lab=1.5)
dev.off()

# Trend, period estimation with estimated period
X=calx(theta1,N)
est.hat <- wsx.bibasis2 (K1,K2,K3,fdPs,fdPt, m,N, ymat,mbasismat, basismat,X)

Trend.hat.Mat <- est.hat$ghat_mat
Period.hat.Mat <- est.hat$mhat_mat
#estgcv=est.hat$gcv

yhat<-(est.hat$ghat+est.hat$mhat)
reshat<-(ymat-est.hat$ghat-est.hat$mhat)

yhat.Mat<-matrix(yhat,m,N)
reshat.Mat<-matrix(reshat,m,N)

Mm <- mesh (s,t)
mx <- Mm$x; my <- Mm$y
TAd=as.Date(TA)
#pdf("hknoxhat.pdf")
#windows(width=14,height=10)
surf3D ( mx,my, yhat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         #zlim=c(zg1,zg2),
         zlab='CO',   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7)


dev.off()

pdf("hknoxres.pdf")
#windows(width=14,height=10)
surf3D ( mx,my, reshat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         #zlim=c(zg1,zg2),
         zlab='CO',   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7)

dev.off()

pdf("ghknoxrest.pdf")
# #windows(width=14,height=10)
plot(TAd,colMeans(reshat.Mat),'l',lty=1,col=1,ylab=expression(paste(epsilon)),xlab='t',cex.axis=1.5,cex.lab=1.5)
dev.off()

pdf("acfnox1.pdf")
#windows(width=10,height=14)
par(mfrow=c(2,1))
acf(colMeans(reshat.Mat),main='ACF')
pacf(colMeans(reshat.Mat),main='PACF')
dev.off()

Box.test(colMeans(reshat.Mat),lag=12)
rem=arima(colMeans(reshat.Mat),order=c(1,0,0))
rem
rem2=arima(colMeans(reshat.Mat),order=c(0,0,1))
rem2

pdf("9ghknoxress.pdf")
#windows(width=14,height=10)
plot(s*24,rowMeans(reshat.Mat),'l',lty=1,col=1,ylab=expression(paste(epsilon)),xlab='hour',cex.axis=1.5,cex.lab=1.5)
dev.off()


Mm <- mesh (s,t_m)
mx <- Mm$x; my <- Mm$y

zg1=min(Trend.hat.Mat)
zg2=max(Trend.hat.Mat)

#pdf("9ghknox.pdf")
#windows(width=14,height=10)
surf3D ( mx,my, Trend.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2),    zlab='g',   theta = 50, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         scale =TRUE, expand =.7,cex.axis=1.2,cex.lab=1.2)

dev.off()
datebreaks <- seq(as.Date("2019-09-01"), as.Date("2021-09-30"),
                  
                  by = "3 month")
# tad=scale_x_date(name = 'time',
#   date_breaks = "3 month",
#   limits = as.Date(c('2019-09-01','2021-09-30')))
TAd=as.Date(TA)


pdf("9ghknoxline.pdf")
#windows(width=7,height=7)
plot(TAd,colMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1,zg2),lty=1,col=1,ylab=expression(paste('g(',s[0],',t/T)')),xlab='t(day)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(TAd,Trend.hat.Mat[3,],lty=2,col=2,lwd=2)
lines(TAd,Trend.hat.Mat[6,],lty=3,col=3,lwd=2)
lines(TAd,Trend.hat.Mat[9,],lty=4,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(TAd,Trend.hat.Mat[12,],lty=5,col=5,lwd=2)
lines(TAd,Trend.hat.Mat[15,],lty=6,col=6,lwd=2)
lines(TAd,Trend.hat.Mat[18,],lty=7,col=7,lwd=2)
lines(TAd,Trend.hat.Mat[21,],lty=8,col=8,lwd=2)
lines(TAd,Trend.hat.Mat[24,],lty=9,col=9,lwd=2)
axis(1,datebreaks,labels=datebreaks)
# lines(t,Trend.hat.Mat[20,],lty=10,col=10)
# lines(t,Trend.hat.Mat[22,],lty=11,col=11)
# lines(t,Trend.hat.Mat[24,],lty=13,col=13)
legend(x=datebreaks[3],y=150,
       c(expression(paste(bar(g),'(t/T)')),expression(paste(s[0],'=3:00')),
         expression(paste(s[0],'=6:00')),
         expression(paste(s[0],'=9:00')),
         expression(paste(s[0],'=12:00')),
         expression(paste(s[0],'=15:00')),
         expression(paste(s[0],'=18:00')),
         expression(paste(s[0],'=21:00')),
         expression(paste(s[0],'=24:00'))),
       lwd=c(2,2,2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7,8,9),
       col=c(1,2,3,4,5,6,7,8,9),cex=1, bty = "n")
dev.off()


pdf("9ghknoxsweek.pdf")
#windows(width=7,height=7)
plot(s*24,rowMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1+13,zg2-4),lty=1,col=1,ylab= expression(paste('g(s,',t[0],'/T)')),xlab='s(hour)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(s*24,Trend.hat.Mat[,1],lty=4,col=2,lwd=3)
lines(s*24,Trend.hat.Mat[,93],lty=13,col=3,lwd=2)
lines(s*24,Trend.hat.Mat[,185],lty=7,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(s*24,Trend.hat.Mat[,277],lty=5,col=5,lwd=2)
lines(s*24,Trend.hat.Mat[,460],lty=6,col=8,lwd=2)
lines(s*24,Trend.hat.Mat[,552],lty=8,col=7,lwd=2)
#lines(s*24,Trend.hat.Mat[,640],lty=8,col=8)
lines(s*24,Trend.hat.Mat[,735],lty=3,col=6,lwd=3)
axis(1,seq(1:24))
legend(x=12,y=65,c(expression(paste(bar(g),'(s)')),expression(paste(t[0],'=1(19-0901-Sun)')),
                   expression(paste(t[0],'=93(19-1202-Mon)')),
                   expression(paste(t[0],'=185(20-0303-Tue)')),
                   expression(paste(t[0],'=277(20-0603-Wed)')),
                   expression(paste(t[0],'=460(20-1203-Thu)')),
                   expression(paste(t[0],'=552(21-0305-Fri)')),
                   expression(paste(t[0],'=735(21-0904-Sat)'))),
       lwd=c(2,3,2,2,2,2,2,3),lty=c(1,4,13,7,5,6,8,3),
       col=c(1,2,3,4,5,8,7,6),cex=1, bty = "n")

dev.off()

Mm1 <- mesh (s,t_m[1:theta1])
mx <- Mm1$x; my <- Mm1$y
zm1=min(Period.hat.Mat[,1:theta1])
zm2=max(Period.hat.Mat[,1:theta1])

#pdf("9mhknox60.pdf")
#windows(width=14,height=10)

surf3D ( mx,my, Period.hat.Mat[,1:theta1], colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  zlab='m',theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7,cex.axis=1.2,cex.lab=1.2)
dev.off()


pdf("9mhknoxline.pdf")
#windows(width=7,height=7)
plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'b',ylim=c(zm1+22,zm2),lty=1,col=1,ylab=expression(paste('m(',s[0],',t)')),
     xlab='t(day)',cex.axis=1.2,cex.lab=1.2,lwd=2,xaxt="n",pch=16)
lines(t_m[1:theta1],Period.hat.Mat[3,1:theta1],lty=3,col=2,lwd=2,'b',pch=9)
lines(t_m[1:theta1],Period.hat.Mat[6,1:theta1],lty=3,col=3,lwd=2,'b',pch=3)
lines(t_m[1:theta1],Period.hat.Mat[9,1:theta1],lty=3,col=4,lwd=2,'b',pch=4)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(t_m[1:theta1],Period.hat.Mat[12,1:theta1],lty=3,col=5,lwd=2,'b',pch=5)
lines(t_m[1:theta1],Period.hat.Mat[15,1:theta1],lty=3,col=6,lwd=2,'b',pch=6)
lines(t_m[1:theta1],Period.hat.Mat[18,1:theta1],lty=3,col=7,lwd=2,'b',pch=11)
lines(t_m[1:theta1],Period.hat.Mat[21,1:theta1],lty=3,col=8,lwd=2,'b',pch=8)
lines(t_m[1:theta1],Period.hat.Mat[24,1:theta1],lty=3,col=9,lwd=2,'b',pch=2)
# lines(t_m[1:theta1],Period.hat.Mat[20,1:theta1],lty=10,col=10)
# lines(t_m[1:theta1],Period.hat.Mat[20,1:theta1],lty=11,col=11)
# lines(t_m[1:theta1],Period.hat.Mat[24,1:theta1],lty=13,col=13)
axis(1, at=seq(1, 7, 1),labels=c('Sun','Mon','Tue','Wed','Thu','Fri','Sat'),cex.axis=1,cex.lab=1,lwd=1)
legend(x=0.8,y=157,c(expression(paste(bar(m),'(t)')),expression(paste(s[0],'=3:00')),
                     expression(paste(s[0],'=6:00')),
                     expression(paste(s[0],'=9:00')),
                     expression(paste(s[0],'=12:00')),
                     expression(paste(s[0],'=15:00')),
                     expression(paste(s[0],'=18:00')),
                     expression(paste(s[0],'=21:00')),
                     expression(paste(s[0],'=24:00'))),
       lwd=c(2,2,2,2,2,2,2,2,2),lty=c(3,3,3,3,3,3,3,3,3), pch=c(16,9,3,4,5,6,11,8,2),
       col=c(1,2,3,4,5,6,7,8,9),cex=1, bty = "n")
dev.off()



pdf("9mhknoxs.pdf")
#windows(width=7,height=7)
plot(s*24,rowMeans(Period.hat.Mat),'l',xaxt="n",ylim=c(zm1,zm2),lty=1,col=1,
     ylab= expression(paste('m(s,',t[0],')')),xlab='s(hour)',cex.axis=1.2,cex.lab=1.2,lwd=2,xaxt="n")
lines(s*24,Period.hat.Mat[,1],lty=4,col=2,lwd=3)
lines(s*24,Period.hat.Mat[,2],lty=13,col=3,lwd=2)
lines(s*24,Period.hat.Mat[,3],lty=7,col=4,lwd=2)
lines(s*24,Period.hat.Mat[,4],lty=5,col=5,lwd=2)
lines(s*24,Period.hat.Mat[,5],lty=6,col=8,lwd=2)
lines(s*24,Period.hat.Mat[,6],lty=8,col=7,lwd=2)
lines(s*24,Period.hat.Mat[,7],lty=3,col=6,lwd=3)
axis(1,seq(1:24))
legend(x=0.1,y=157,c(expression(paste(bar(m),'(s)')),expression(paste(t[0],'=1+7k(Sun)')),
                     expression(paste(t[0],'=2+7k(Mon)')),expression(paste(t[0],'=3+7k(Tue)')),
                     expression(paste(t[0],'=4+7k(Wed)')),expression(paste(t[0],'=5+7k(Thu)')),
                     expression(paste(t[0],'=6+7k(Fri)')),expression(paste(t[0],'=7+7k(Sat)'))),
       lwd=c(2,3,2,2,2,2,2,3),lty=c(1,4,13,7,5,6,8,3),
       col=c(1,2,3,4,5,8,7,6),cex=1, bty = "n")

dev.off()


m7=c(Period.hat.Mat[,1],Period.hat.Mat[,2],Period.hat.Mat[,3],
     Period.hat.Mat[,4],Period.hat.Mat[,5],Period.hat.Mat[,6]
     ,Period.hat.Mat[,7])
s7=seq(1,168,1)
pdf("9mhknoxsweek.pdf")
#windows(width=7,height=7)
plot(s7,m7,'l',
     xaxt="n",ylim=c(zm1,zm2),lty=1,col=1,ylab='m',
     xlab='time',cex.axis=1.5,cex.lab=1.5,lwd=2,ann = F)
abline(v=24,lty=4)
abline(v=48,lty=4)
abline(v=72,lty=4)
abline(v=96,lty=4)
abline(v=120,lty=4)
abline(v=144,lty=4)
abline(v=168,lty=4)
axis(1,  tck=0,at=seq(12, 168, 24),labels=c('Sun','Mon','Tue','Wed','Thu','Fri','Sat'),cex.axis=1.2,cex.lab=1.2,lwd=1)
title( ylab = 'm',cex.axis=1.5,cex.lab=1.5)


dev.off()


Mm1 <- mesh (s,t_m)
mx <- Mm1$x; my <- Mm1$y
#pdf("9pmhknox.pdf")
#windows(width=14,height=10)
surf3D ( mx,my, Period.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  zlab='m',  theta = 45, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         #main = TeX ('\\hat{m}(s,t)'),
         scale =TRUE , expand =.7 )
dev.off()

Mm <- mesh (s,t_m)
mx <- Mm$x; my <- Mm$y
zg1=min(min(YY),min(YY))
zg2=max(max(YY),max(YY))

#pdf("hknox.pdf")
#windows(width=14,height=10)
surf3D ( mx,my, YY, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2), zlab='NOX',   theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7 ,cex.axis=1.2,cex.lab=1.2)


dev.off()

