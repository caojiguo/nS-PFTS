# this function differs with the original one in the penalties:
# we consider 
# lambdat*kronecker(penmatt,sinprod) instead of 
# lambdas*kronecker(Imat,penmats) for s direction
# and, we consider 
# lambdat*kronecker(penmatt,sinprod) instead of
# lambdat*kronecker(penmatt,Imat) for t direction 
# the modification is according to marginal penalization: Wood 2006 (see notes)
wsx.bibasis <- function (ymat,mbasismat, basismat,X){

  basismatx  = kronecker(X,mbasismat)
  basismata = cbind(basismat,basismatx)
  
  Bmat  = crossprod(basismata,basismata)
  
  Dmat = crossprod(basismata,ymat)
  Bmatinv=ginv(Bmat)
  coef =Bmatinv%*% Dmat 
  
  yhat = basismata %*% coef
  SSE  = sum((ymat - yhat)^2)
  
  #  class(smoothlist) = "bifdSmooth"
  return(SSE)
  
}

wsx.bibasis2 <- function (K1,K2,K3,fdPars,fdPart, m,N, ymat,mbasismat, basismat,X,fdnames=NULL){
  
  if (is.null(fdnames)) {
    fdnames      = vector("list", 3)
    fdnames[[1]] = "argument s"
    fdnames[[2]] = "argument t"
    fdnames[[3]] = "function"
  }
  xdim=dim(X)
  fdobjs  = fdPars$fd
  sbasis  = fdobjs$basis
  snbasis = sbasis$nbasis - length(sbasis$dropind)


  
  #fdPart  = fdParcheck(fdPart)
  fdobjt  = fdPart$fd
  tbasis  = fdobjt$basis
  tnbasis = tbasis$nbasis - length(tbasis$dropind)

  basismatx  = kronecker(X,mbasismat)
  basismata = cbind(basismat,basismatx)
  
  Bmat  = crossprod(basismata,basismata)
  
  Dmat = crossprod(basismata,ymat)
  Bmatinv=ginv(Bmat)
  coef =Bmatinv%*% Dmat 
 
  coefmatm = matrix(coef[(K1*K2+1):length(coef)], K3, xdim[2])
  coefmat = matrix(coef[1:(K1*K2)], K1, K2)
  bifdobj = bifd(coefmat, sbasis, tbasis, fdnames)
  yhat = basismata %*% coef
  SSE  = sum((ymat - yhat)^2)
  ghat= basismat%*% coef[1:(K1*K2)]
  ghat_mat=matrix(ghat,m,N)
  # 
  mhat=basismatx%*% coef[(K1*K2+1):length(coef)]
  mhat_mat=matrix(mhat,m,N)
  #  compute  GCV index  
  
  BiB0 = Bmatinv %*% Bmat
  
  df = sum(diag(BiB0))
  Nt = (m*N)
  if (df < Nt) {
    gcv = (SSE/Nt)/((Nt - df)/Nt)^2
  } else {
    gcv = NA
  }
  
  smoothlist = list(df=df,     bifdobj=bifdobj,      gcv=gcv,
                    SSE=SSE,  coff=coef,  coefg = coefmat,coefm=coefmatm,
                    ghat_mat=ghat_mat, mhat_mat=mhat_mat,ghat=ghat,mhat=mhat)
  
  #  class(smoothlist) = "bifdSmooth"
  return(smoothlist)
  
}

#que no regerse penmats