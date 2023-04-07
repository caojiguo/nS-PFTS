require ( plot3D )
require (fda )
require ( latex2exp )
require ( pracma )
require ( mgcv )

setwd("E:/paper/functional period/code")

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

mdata <- read.csv("global1850.csv")

ymat<-as.matrix(mdata[,1])
#ymat<-as.matrix(mdata[,8])
m<-12
N<-length(ymat)/m
YY<-matrix(ymat,m,N)

TAd<-seq (1850,2021, by=1)

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
Thm_0=100
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
  Q[k]=Q2[k]+log(N*m)*(thetai)*s1*(N*m)^0.4/6
  k=k+1
}
theta1 = pseq[which.min(Q)]
# pdf("Qgltem.pdf")
#windows(width=7,height=7)
plot(pseq,Q,"l",ylab='Q',xlab=expression(paste(theta)),cex.axis=1.5,cex.lab=1.5,lwd=2)
axis(1,theta1,labels=theta1,cex.axis=1.5,cex.lab=1.5)
#
dev.off()

# Trend, period estimation with estimated period
X=calx(theta1,N)
est.hat <- wsx.bibasis2 (K1,K2,K3,fdPs,fdPt, m,N, ymat,mbasismat, basismat,X)

Trend.hat.Mat <- est.hat$ghat_mat
Period.hat.Mat <- est.hat$mhat_mat
#estgcv=est.hat$gcv



Mm <- mesh (s,TAd)
mx <- Mm$x; my <- Mm$y

zg1=min(Trend.hat.Mat)
zg2=max(Trend.hat.Mat)

#pdf("gglobal.pdf")
#windows(width=7,height=7)
surf3D ( mx,my, Trend.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",cex.axis=1.2,cex.lab=1.2,
         zlim=c(zg1,zg2),    zlab='g',   theta = 50, phi =30, facets = FALSE,bty="g", ticktype = "detailed",
         scale =TRUE, expand =.7 
         #,plot=FALSE,main = TeX ('$\\hat{g}(s,t)')
)

dev.off()

pdf("gtemline12.pdf")
#windows(width=7,height=7)
plot(TAd,colMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1,zg2),lty=1,col=1,
     ylab=expression(paste('g(',s[0],',t/T)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(TAd,Trend.hat.Mat[1,],lty=2,col=13,lwd=2)
lines(TAd,Trend.hat.Mat[2,],lty=3,col=3,lwd=2)
lines(TAd,Trend.hat.Mat[3,],lty=4,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(TAd,Trend.hat.Mat[4,],lty=5,col=5,lwd=2)
lines(TAd,Trend.hat.Mat[5,],lty=6,col=6,lwd=2)
lines(TAd,Trend.hat.Mat[6,],lty=12,col=12,lwd=2)
lines(TAd,Trend.hat.Mat[7,],lty=8,col=8,lwd=2)
lines(TAd,Trend.hat.Mat[8,],lty=9,col=9,lwd=2)
lines(TAd,Trend.hat.Mat[9,],lty=10,col=10,lwd=2)
lines(TAd,Trend.hat.Mat[10,],lty=11,col=11,lwd=2)
lines(TAd,Trend.hat.Mat[11,],lty=7,col=7,lwd=3)
lines(TAd,Trend.hat.Mat[12,],lty=13,col=2,lwd=3)
axis(1,seq(1850,2019,by=13),labels=seq(1850,2019,by=13))
# lines(t,Trend.hat.Mat[20,],lty=10,col=10)
# lines(t,Trend.hat.Mat[22,],lty=11,col=11)
# lines(t,Trend.hat.Mat[24,],lty=13,col=13)
legend(x=1863,y=1,
       c(expression(paste(bar(g),'(t/T)')),expression(paste(s[0],'=1(Jan)')),
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
       lwd=c(2,2,2,2,2,2,2,2,2,2,2,3,3),lty=c(1,2,3,4,5,6,12,8,9,10,11,7,13),
       col=c(1,13,3,4,5,6,12,8,9,10,11,7,2),cex=1, bty = "n")
dev.off()




pdf("gtemsline.pdf")
#windows(width=7,height=7)
plot(s*12,rowMeans(Trend.hat.Mat),'l',xaxt="n",ylim=c(zg1,zg2),
     #     ylim=c(min(Trend.hat.Mat[,1],rowMeans(Trend.hat.Mat)),max(Trend.hat.Mat[,1],rowMeans(Trend.hat.Mat))),
     lty=1,col=1,ylab=expression(paste('g(s,',t[0],'/T)')),xlab='s(month)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(s*12,Trend.hat.Mat[,1],lty=2,col=2,lwd=2)
lines(s*12,Trend.hat.Mat[,41],lty=3,col=3,lwd=2)
lines(s*12,Trend.hat.Mat[,81],lty=4,col=4,lwd=2)
lines(s*12,Trend.hat.Mat[,121],lty=5,col=5,lwd=2)
lines(s*12,Trend.hat.Mat[,151],lty=6,col=6,lwd=2)
lines(s*12,Trend.hat.Mat[,161],lty=7,col=7,lwd=2)
lines(s*12,Trend.hat.Mat[,171],lty=8,col=8,lwd=2)
axis(1,seq(1:12))
legend(x=9.9,y=0.49,
       c(expression(paste(bar(g),'(s)')),expression(paste(t[0],'=1850')),
         expression(paste(t[0],'=1890')),
         expression(paste(t[0],'=1930')),
         expression(paste(t[0],'=1970')),
         expression(paste(t[0],'=2000')),
         expression(paste(t[0],'=2010')),
         expression(paste(t[0],'=2020'))),
       lwd=c(2,2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7,8),
       col=c(1,2,3,4,5,6,7,8),cex=1, bty = "n")
dev.off()


Mm <- mesh (s,TAd)
mx <- Mm$x; my <- Mm$y
zg1=min(min(YY),min(YY))
zg2=max(max(YY),max(YY))

pdf("global.pdf",width=8,height=7)
#windows(width=8,height=7)
#par(mgp=c(3,2,0))
surf3D ( mx,my, YY, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zg1,zg2), zlab='Tem',   theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
     scale =TRUE , expand =.7,cex.axis=1.5,cex.lab=1.5)

dev.off()

Mm1 <- mesh (s,t_m[1:theta1])
mx <- Mm1$x; my <- Mm1$y
zm1=min(Period.hat.Mat[,1:theta1])
zm2=max(Period.hat.Mat[,1:theta1])

pdf("temm.pdf")
#windows(width=7,height=7)

surf3D ( mx,my, Period.hat.Mat[,1:theta1], colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
         zlim=c(zm1,zm2),  zlab='m',theta = 50, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
         scale =TRUE , expand =.7 ,cex.axis=1.2,cex.lab=1.2)
dev.off()



# pdf("mtemline16.pdf")
# #windows(width=7,height=7)
# plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'l',ylim=c(zm1,zm2),
#      lty=1,col=1,ylab=expression(paste('m(',s[0],',t)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[1,1:theta1],lty=2,col=2,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[2,1:theta1],lty=3,col=3,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[3,1:theta1],lty=4,col=4,lwd=2)
# #lines(t,Trend.hat.Mat[13,],lty=5,col=5)
# lines(t_m[1:theta1],Period.hat.Mat[4,1:theta1],lty=5,col=5,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[5,1:theta1],lty=6,col=6,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[6,1:theta1],lty=7,col=7,lwd=2)
# 
# legend(x=14,y=0.38,c(expression(paste(bar(m),'(t)')),expression(paste(s[0],'=1(Jan)')),
#                      expression(paste(s[0],'=2(Feb)')),
#                      expression(paste(s[0],'=3(Mar)')),
#                      expression(paste(s[0],'=4(Apr)')),
#                      expression(paste(s[0],'=5(May)')),
#                      expression(paste(s[0],'=6(June)'))),
#        lwd=c(2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7),
#        col=c(1,2,3,4,5,6,7),cex=1, bty = "n")
# dev.off()
# 
# pdf("mtemline712.pdf")
# #windows(width=7,height=7)
# plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'l',ylim=c(zm1,zm2),
#      lty=1,col=1,ylab=expression(paste('m(',s[0],',t)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[7,1:theta1],lty=2,col=2,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[8,1:theta1],lty=3,col=3,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[9,1:theta1],lty=4,col=4,lwd=2)
# #lines(t,Trend.hat.Mat[13,],lty=5,col=5)
# lines(t_m[1:theta1],Period.hat.Mat[10,1:theta1],lty=5,col=5,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[11,1:theta1],lty=6,col=6,lwd=2)
# lines(t_m[1:theta1],Period.hat.Mat[12,1:theta1],lty=7,col=7,lwd=2)
# 
# legend(x=14,y=0.38,c(expression(paste(bar(m),'(t)')),expression(paste(s[0],'=7(July)')),
#        expression(paste(s[0],'=8(Aug)')), 
#        expression(paste(s[0],'=9(Sept)')),
#        expression(paste(s[0],'=10(Oct)')),
#        expression(paste(s[0],'=11(Nov)')),
#        expression(paste(s[0],'=12(Dec)'))),
#        lwd=c(2,2,2,2,2,2,2),lty=c(1,2,3,4,5,6,7),
#        col=c(1,2,3,4,5,6,7),cex=1, bty = "n")
# dev.off()



pdf("mtemline12.pdf")
#windows(width=7,height=7)
plot(t_m[1:theta1],colMeans(Period.hat.Mat)[1:theta1],'l',ylim=c(zm1,zm2),lty=1,col=1,
     ylab=expression(paste('m(',s[0],',t)')),xlab='t(year)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[1,1:theta1],lty=7,col=12,lwd=3)
lines(t_m[1:theta1],Period.hat.Mat[2,1:theta1],lty=3,col=3,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[3,1:theta1],lty=4,col=4,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(t_m[1:theta1],Period.hat.Mat[4,1:theta1],lty=5,col=5,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[5,1:theta1],lty=6,col=6,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[6,1:theta1],lty=2,col=7,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[7,1:theta1],lty=8,col=8,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[8,1:theta1],lty=9,col=9,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[9,1:theta1],lty=10,col=10,lwd=2)
#lines(t,Trend.hat.Mat[13,],lty=5,col=5)
lines(t_m[1:theta1],Period.hat.Mat[10,1:theta1],lty=11,col=11,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[11,1:theta1],lty=12,col=13,lwd=2)
lines(t_m[1:theta1],Period.hat.Mat[12,1:theta1],lty=13,col=2,lwd=3)
legend(x=19,y=0.39, c(expression(paste(bar(m),'(t)')),expression(paste(s[0],'=1')),
                      expression(paste(s[0],'=2')),
                      expression(paste(s[0],'=3')),
                      expression(paste(s[0],'=4')),
                      expression(paste(s[0],'=5')),
                      expression(paste(s[0],'=6')),
                      expression(paste(s[0],'=7')),
                      expression(paste(s[0],'=8')), 
                      expression(paste(s[0],'=9')),
                      expression(paste(s[0],'=10')),
                      expression(paste(s[0],'=11')),
                      expression(paste(s[0],'=12'))),
       lwd=c(2,3,2,2,2,2,2,2,2,2,2,2,3),lty=c(1,7,3,4,5,6,2,8,9,10,11,12,13),
       col=c(1,12,3,4,5,6,7,8,9,10,11,13,2),cex=0.9, bty = "n")
dev.off()



pdf("mtems2.pdf")
#windows(width=7,height=7)
plot(s*12,rowMeans(Period.hat.Mat),'l',xaxt="n",ylim=c(zm1,zm2),lty=1,col=1, 
     ylab= expression(paste('m(s,',t[0],')')),xlab='s(month)',cex.axis=1.2,cex.lab=1.2,lwd=2)
lines(s*12,Period.hat.Mat[,1],lty=2,col=2,lwd=2)
lines(s*12,Period.hat.Mat[,13],lty=9,col=9,lwd=2)
lines(s*12,Period.hat.Mat[,23],lty=4,col=4,lwd=2)
lines(s*12,Period.hat.Mat[,29],lty=5,col=5,lwd=2)
lines(s*12,Period.hat.Mat[,33],lty=6,col=6,lwd=2)
lines(s*12,Period.hat.Mat[,35],lty=7,col=7,lwd=2)
lines(s*12,Period.hat.Mat[,41],lty=8,col=8,lwd=2)
# lines(s*12,Period.hat.Mat[,7],lty=8,col=8)
axis(1,seq(1:12))
legend(x=8,y=0.37,c(expression(paste(bar(m),'(s)')),
                    expression(paste(t[0],'=1850+46k')),
                    expression(paste(t[0],'=1862+46k')),expression(paste(t[0],'=1872+46k')),
                    expression(paste(t[0],'=1878+46k')),expression(paste(t[0],'=1882+46k')),
                    expression(paste(t[0],'=1884+46k')),expression(paste(t[0],'=1890+46k'))),
       # c("mean","1850,1896,1942,1988","1862,1908,1954,2000",
       #              "1872,1918,1964,2010",
       #              "1878,1924,1970,2016","1882,1928,1974,2020",
       #              "1884,1930,1976","1890,1936,1982"),
       lwd=c(2,2,2,2,2,2,2,2),lty=c(1,2,9,4,5,6,7,8),
       col=c(1,2,9,4,5,6,7,8),cex=1, bty = "n")

dev.off()

# 
# Mm1 <- mesh (s,t_m)
# mx <- Mm1$x; my <- Mm1$y
# pdf("ptemm.pdf")
# #windows(width=7,height=7)
# surf3D ( mx,my, Period.hat.Mat, colkey = TRUE,cex.main =2, xlab ="s", ylab ="t",
#          zlim=c(zm1,zm2),  zlab='m',  theta = 60, phi =30, facets = FALSE ,bty="g", ticktype = "detailed",
#          #main = TeX ('\\hat{m}(s,t)'),
#          scale =TRUE , expand =.7 )
# dev.off()


hkco_gcv<-function(N,m,ymat,theta0,K11,K22,K33,c_a){
  gcv_bic=matrix(0,(K11-3)*(K22-3)*(K33-4),2)
  kseq=matrix(0,(K11-3)*(K22-3)*(K33-4),3)
  s <- seq (1,m, by=1)/m
  t <- seq (1,N, by=1)/N
  t_m <- seq (1,N, by=1)
  # create basis for each coordinates
  
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

theta0=theta1=46
K11=10
K22=10
K33=10
c_a=1
s_time=proc.time()
resultgcv<-hkco_gcv(N,m,ymat,theta0,K11,K22,K33,c_a)
e_time=proc.time()

c(e_time-s_time)

c(which.min(resultgcv$gcv_bic[,1]),which.min(resultgcv$gcv_bic[,2]))

resultgcv$kseq[which.min(resultgcv$gcv_bic[,1]),]

resultgcv$kseq[which.min(resultgcv$gcv_bic[,2]),]

bicc=as.matrix(resultgcv$gcv_bic[,2],1)
gcvv=as.matrix(resultgcv$gcv_bic[,1],1)
ind_k=seq(1,length(bicc),by=1)
windows(width=7,height=7)
plot(ind_k,bicc,'l')

windows(width=7,height=7)
plot(ind_k,gcvv,'l')

resultgcv$kseq[which.min(resultgcv$gcv_bic[,1]),]

resultgcv$kseq[which.min(resultgcv$gcv_bic[,2]),]
