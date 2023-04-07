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

mdata <- read.csv("sunspots.csv")
TA=as.matrix(mdata[1:2820,1])
TAd<-seq (1749,1983, by=1)
ymat<-as.matrix(mdata[,2])
#ymat<-as.matrix(mdata[,8])
m<-12
N<-length(ymat)/m
YY<-matrix(ymat,m,N)

s <- seq (1,m, by=1)/m
t <- seq (1,N, by=1)/N
t_m <- seq (1,N, by=1)
K1=6
K3=6
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
  Q[k]=Q2[k]+c_a*log(N*m)*(thetai)*s1*(N*m)^(1/6)
  k=k+1
}
theta1 = pseq[which.min(Q)]
theta1
pdf("Qsun.pdf")
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

#pdf("sunhat.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, yhat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         #zlim=c(zg1,zg2),
         zlab='sunspots',   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7)


dev.off()

pdf("sunres.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, reshat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         #zlim=c(zg1,zg2),
         zlab='sunspots',   theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7)

dev.off()

pdf("sunrest.pdf")
# #windows(width=7,height=7)
plot(TAd,colMeans(reshat.Mat),'l',lty=1,col=1,ylab=expression(paste(epsilon)),xlab='t',cex.axis=1.5,cex.lab=1.5)

dev.off()

pdf("acfsun1.pdf")
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

pdf("sunress.pdf")
#windows(width=7,height=7)
plot(s*12,rowMeans(reshat.Mat),'l',lty=1,col=1,ylab=expression(paste(epsilon)),xlab='hour',cex.axis=1.5,cex.lab=1.5)
dev.off()


Mm <- mesh (s,TAd)
mx <- Mm$x; my <- Mm$y

zg1=min(Trend.hat.Mat)
zg2=max(Trend.hat.Mat)

#pdf("suntrend.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, Trend.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",cex.axis=1.2,cex.lab=1.2,
         zlim=c(zg1,zg2),    zlab='g',   theta = 50, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         scale =TRUE, expand =.7 
         #,plot=FALSE,main = TeX ('$\\hat{g}(s,t)')
)

dev.off()

pdf("gsunline12.pdf")
#windows(width=7,height=7)
plot(TAd,colMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1,zg2),lty=1,col=1,
     ylab=expression(paste('g(',s[0],',t/T)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(TAd,Trend.hat.Mat[1,],lty=2,col=2,lwd=2)
lines(TAd,Trend.hat.Mat[2,],lty=3,col=3,lwd=2)
lines(TAd,Trend.hat.Mat[3,],lty=4,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(TAd,Trend.hat.Mat[4,],lty=5,col=5,lwd=2)
lines(TAd,Trend.hat.Mat[5,],lty=6,col=6,lwd=2)
lines(TAd,Trend.hat.Mat[6,],lty=7,col=7,lwd=2)
lines(TAd,Trend.hat.Mat[7,],lty=8,col=8,lwd=2)
lines(TAd,Trend.hat.Mat[8,],lty=9,col=9,lwd=2)
lines(TAd,Trend.hat.Mat[9,],lty=10,col=10,lwd=2)
lines(TAd,Trend.hat.Mat[10,],lty=11,col=11,lwd=2)
lines(TAd,Trend.hat.Mat[11,],lty=12,col=12,lwd=2)
lines(TAd,Trend.hat.Mat[12,],lty=13,col=13,lwd=2)
axis(1,seq(1749,1983,by=6),labels=seq(1749,1983,by=6))
#axis(1,seq(1749,1983,by=1),labels=seq(1749,1983,by=1))
# lines(t,Trend.hat.Mat[20,],lty=10,col=10)
# lines(t,Trend.hat.Mat[22,],lty=11,col=11)
# lines(t,Trend.hat.Mat[24,],lty=13,col=13)
legend(x=1865,y=73, c(expression(paste(bar(g),'(t/T)')),expression(paste(s[0],'=1(Jan)')),
                      expression(paste(s[0],'=2(Feb)')),
                      expression(paste(s[0],'=3(Mar)')),
                      expression(paste(s[0],'=4(Apr)')),
                      expression(paste(s[0],'=5(May)')),
                      expression(paste(s[0],'=6(June)')),
                      expression(paste(s[0],'=7(July)')),
                      expression(paste(s[0],'=8(Aug)')), 
                      expression(paste(s[0],'=9(Sept)')),
                      expression(paste(s[0],'=10(Oct)')),
                      expression(paste(s[0],'=11(Nov)')),
                      expression(paste(s[0],'=12(Dec)'))),
       lwd=c(2,2,2,2,2,2,2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7,8,9,10,11,12,13),
       col=c(1,2,3,4,5,6,7,8,9,10,11,12,13),cex=1, bty = "n")
dev.off()



pdf("gsunsline.pdf")
#windows(width=7,height=7)
plot(s*12,rowMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1,zg2),
     #     ylim=c(min(Trend.hat.Mat[,1],rowMeans(Trend.hat.Mat)),max(Trend.hat.Mat[,1],rowMeans(Trend.hat.Mat))),
     lty=1,col=1,ylab=expression(paste('g(s,',t[0],'/T)')),xlab='s(month)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(s*12,Trend.hat.Mat[,1],lty=2,col=2,lwd=2)
lines(s*12,Trend.hat.Mat[,35],lty=3,col=3,lwd=2)
lines(s*12,Trend.hat.Mat[,69],lty=4,col=4,lwd=2)
lines(s*12,Trend.hat.Mat[,103],lty=5,col=5,lwd=2)
lines(s*12,Trend.hat.Mat[,137],lty=6,col=6,lwd=2)
lines(s*12,Trend.hat.Mat[,171],lty=7,col=7,lwd=2)
lines(s*12,Trend.hat.Mat[,205],lty=8,col=8,lwd=2)
lines(s*12,Trend.hat.Mat[,235],lty=9,col=9,lwd=2)
axis(1,seq(1:12))
legend(x=10.2,y=73.2,c(expression(paste(bar(g),'(s)')),expression(paste(t[0],'=1749')),
                   expression(paste(t[0],'=1783')),
                   expression(paste(t[0],'=1817')),
                   expression(paste(t[0],'=1851')),
                   expression(paste(t[0],'=1885')),
                   expression(paste(t[0],'=1919')),
                   expression(paste(t[0],'=1953')),
                   expression(paste(t[0],'=1983'))),
       lwd=c(2,2,2,2,2,2,2,2,2),
       # c("mean","1749","1783","1817","1851","1885","1919","1953","1983"),
       # lwd=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5),
       lty=c(1,2,3,4,5,6,7,8,9),
       col=c(1,2,3,4,5,6,7,8,9),cex=0.9, bty = "n")
dev.off()


Mm <- mesh (s,TAd)
mx <- Mm$x; my <- Mm$y
zg1=min(min(YY),min(YY))
zg2=max(max(YY),max(YY))

#pdf("sun.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, YY, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2), zlab='sunspots',   theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7,cex.axis=1.2,cex.lab=1.2 )


dev.off()

Mm1 <- mesh (s,t_m[1:theta1])
mx <- Mm1$x; my <- Mm1$y
zm1=min(Period.hat.Mat[,1:theta1])
zm2=max(Period.hat.Mat[,1:theta1])

#pdf("sunm.pdf")
#windows(width=7,height=7)

surf3D ( mx,my, Period.hat.Mat[,1:theta1], colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  zlab='m',theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7,cex.axis=1.5,cex.lab=1.5 
         #,plot=FALSE
)

dev.off()



pdf("msunline12.pdf")
#windows(width=7,height=7)
plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'l',ylim=c(zm1,zm2),lty=1,col=1,
     ylab=expression(paste('m(',s[0],',t)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[1,1:theta1],lty=2,col=2,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[2,1:theta1],lty=3,col=3,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[3,1:theta1],lty=4,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(t_m[1:theta1],Period.hat.Mat[4,1:theta1],lty=5,col=5,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[5,1:theta1],lty=6,col=6,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[6,1:theta1],lty=7,col=7,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[7,1:theta1],lty=8,col=8,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[8,1:theta1],lty=9,col=9,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[9,1:theta1],lty=10,col=10,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(t_m[1:theta1],Period.hat.Mat[10,1:theta1],lty=11,col=11,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[11,1:theta1],lty=12,col=12,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[12,1:theta1],lty=13,col=13,lwd=3)
legend(x=6,y=68,c(expression(paste(bar(m),'(t)')),expression(paste(s[0],'=1(Jan)')),
                  expression(paste(s[0],'=2(Feb)')),
                  expression(paste(s[0],'=3(Mar)')),
                  expression(paste(s[0],'=4(Apr)')),
                  expression(paste(s[0],'=5(May)')),
                  expression(paste(s[0],'=6(June)')),
                  expression(paste(s[0],'=7(July)')),
                  expression(paste(s[0],'=8(Aug)')), 
                  expression(paste(s[0],'=9(Sept)')),
                  expression(paste(s[0],'=10(Oct)')),
                  expression(paste(s[0],'=11(Nov)')),
                  expression(paste(s[0],'=12(Dec)'))),
       lwd=c(2,2,2,2,2,2,2,2,2,2,2,2,3),lty=c(1,2,3,4,5,6,7,8,9,10,11,12,13),
       col=c(1,2,3,4,5,6,7,8,9,10,11,12,13),cex=1, bty = "n")
dev.off()
# pdf("mhkcot.pdf")
# #windows(width=7,height=7)
# plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'l',lty=1,col=1,ylab='m',xlab='t',cex.axis=1.5,cex.lab=1.5)
# dev.off()


pdf("msunsy12.pdf")
#windows(width=7,height=7)
plot(s*12,rowMeans(Period.hat.Mat),'l',xaxt="n",ylim=c(zm1,zm2),lty=1,col=1,
     ylab= expression(paste('m(s,',t[0],')')),xlab='s(month)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(s*12,Period.hat.Mat[,1],lty=2,col=2,lwd=2)
lines(s*12,Period.hat.Mat[,2],lty=3,col=3,lwd=2)
lines(s*12,Period.hat.Mat[,3],lty=4,col=4,lwd=2)
lines(s*12,Period.hat.Mat[,4],lty=5,col=5,lwd=2)
lines(s*12,Period.hat.Mat[,5],lty=6,col=6,lwd=2)
lines(s*12,Period.hat.Mat[,6],lty=7,col=7,lwd=2)
lines(s*12,Period.hat.Mat[,7],lty=8,col=8,lwd=2)
 lines(s*12,Period.hat.Mat[,8],lty=9,col=9,lwd=2)
 lines(s*12,Period.hat.Mat[,9],lty=10,col=10,lwd=2)
 lines(s*12,Period.hat.Mat[,10],lty=11,col=11,lwd=2)
 lines(s*12,Period.hat.Mat[,11],lty=12,col=12,lwd=2)
axis(1,seq(1:12))
legend(x=9.1,y=70.5, c(expression(paste(bar(m),'(s)')),
                   expression(paste(t[0],'=1749+11k')),
                   expression(paste(t[0],'=1750+11k')),expression(paste(t[0],'=1751+11k')),
                   expression(paste(t[0],'=1752+11k')),expression(paste(t[0],'=1753+11k')),
                   expression(paste(t[0],'=1754+11k')),expression(paste(t[0],'=1755+11k')),
                   expression(paste(t[0],'=1756+11k')),expression(paste(t[0],'=1757+11k'))),
       # c("mean","1749","1750","1751","1752","1753","1754","1755","1756","1757",
       #            "1758","1759"),
       lwd=c(2,2,2,2,2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7,8,9,10,11,12),
       col=c(1,2,3,4,5,6,7,8,9,10,11,12),cex=1, bty = "n")

dev.off()



Mm1 <- mesh (s,t_m)
mx <- Mm1$x; my <- Mm1$y
pdf("psunm.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, Period.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  zlab='m',  theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         #main = TeX ('\\hat{m}(s,t)'),
         scale =TRUE , expand =.7 )
dev.off()




