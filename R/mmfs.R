################################################################################
##
##   Copyright (C) 2015 - 2016 Alfred Galichon
##
##   This file is part of the R package TraME.
##
##   The R package TraME is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package TraME is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
#
# References:
#
# A. Galichon, S. Weber: "Estimation of Matching Function Equilibria"
#           
#           
# mmfs has n, and m and neededNorm
# M.mmfs <- function(mmfs, ax, by)
#
################################################################################
########################    Default and generic methods ########################
################################################################################
#
margxInv.default <- function(xs, mmfs, theMu0ys)
{
  nbX = length(mmfs$n)
  if (is.null(mmfs$neededNorm))
  {
    coeff = 1
    ubs = mmfs$n
  }
  else
  {
    coeff = 0
    ubs = rep(1e10,nbX)
  }
  #
  if(is.null(xs)){
    xs = 1:nbX
  }
  theaxs = rep(0,length(xs))
  #
  i = 0
  for(x in xs){
    i = i+1
    root_fn <- function(z) (ifelse(coeff,z,0) - mmfs$n[x] + sum(M(mmfs,z,theMu0ys,xs=x)))
    theaxs[i] = uniroot(root_fn, c(0,ubs[x]), tol = 1e-300)$root # Keith: fix tolerence   
  }
  #
  return(theaxs)
}
#
margyInv.default <- function(ys, mmfs, theMux0s)
{
  nbY = length(mmfs$m)
  if (is.null(mmfs$neededNorm))
  {
    coeff = 1
    ubs = mmfs$m
  }
  else
  {
    coeff = 0
    ubs = rep(1e10,nbY)
  }
  #
  if(is.null(ys)){
    ys = 1:nbY
  }
  thebys = rep(0,length(ys))
  #
  j = 0
  for(y in ys){
    j = j+1
    root_fn <- function(z) (ifelse(coeff,z,0) - mmfs$m[y] + sum(M(mmfs,theMux0s,z,ys=y)))
    thebys[j] = uniroot(root_fn, c(0,ubs[y]), tol=1e-300)$root
  }
  #
  return(thebys)
}
#
################################################################################
########################             TU MMfs            ########################
################################################################################
build_TUmmfs <- function(n,m,K,neededNorm)
{
  ret = list(n=n,
             m=m,
             nbX=length(n),
             nbY=length(m),
             neededNorm=neededNorm,
             K = K
  )
  class(ret)="TUmmfs"
  return(ret)
}
#
mmfsTranspose.TUmmfs <- function(mmfs)
{
  ret = list(n=mmfs$m,
             m=mmfs$n,
             nbX = mmfs$nbY,
             nbY = mmfs$nbX,
             neededNorm=normalizationTranspose(mmfs$neededNorm),
             K = t(mmfs$K)
  )
  class(ret)="TUmmfs"
  return(ret)
}
#
M.TUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = mmfs$K[xs,ys]
  term_2 = sqrt(mux0s[xs] %*% t(mu0ys[ys]))
  #
  ret = term_1 * term_2
  #
  return(ret)
}
# 
dmux0s_M.TUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mmfs$K / 2
  term_2 = sqrt( (1/mux0s) %*% t(mu0ys))
  #
  ret = term_1 * term_2
  #
  return(ret)
}
#
dmu0ys_M.TUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mmfs$K / 2
  term_2 = sqrt( mux0s %*% t(1/mu0ys))
  #
  ret = term_1 * term_2
  #
  return(ret)
}
#
dtheta_M.TUmmfs <- function (mmfs,mux0s,mu0ys,dtheta=NULL)
{
  ret <- 0
  if(is.null(dtheta)){
    ret = Diagonal(sqrt(mux0s %*% t(mu0ys)),n=mmfs$nbX*mmfs$nbY)
  }else{
    ret = c(matrix(dtheta,nrow = nbX) * sqrt(mux0s %*% t(mu0ys)) )
  }
  return(ret)
}
#
margxInv.TUmmfs <- function(xs, mmfs, theMu0ys)
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  #
  sqrtBys = sqrt(theMu0ys)
  if(is.null(mmfs$neededNorm)){
    b = (mmfs$K[xs,] %*% sqrtBys)/2
    sqrtAxs = sqrt(mmfs$n[xs]+ b*b) - b
  }else{
    sqrtAxs = mmfs$n / c(mmfs$K[xs,] %*% sqrtBys)
  }
  #
  ret = c(sqrtAxs*sqrtAxs)
  #
  return(ret)
}
#
margyInv.TUmmfs <- function(ys, mmfs, theMux0s)
{
  if (is.null(ys)) {xs = 1:length(mmfs$m)}
  #
  sqrtAxs = sqrt(theMux0s)
  if(is.null(mmfs$neededNorm)){
    b = c(sqrtAxs %*% mmfs$K[,ys])/2
    sqrtBys = sqrt(mmfs$m[ys] + b*b) - b
  }else{
    sqrtBys = mmfs$m / c(sqrtAxs %*% mmfs$K[,ys])
  }
  #
  ret = sqrtBys*sqrtBys
  #
  return(ret)
}

################################################################################
########################            NTU MMfs            ########################
################################################################################
build_NTUmmfs <- function(n,m,A,B,neededNorm)
{
  ret = list(n=n,
             m=m,
             nbX=length(n),
             nbY=length(m),
             neededNorm=neededNorm,
             A = A,
             B = B)
  class(ret)="NTUmmfs"
  return(ret)
}
#
mmfsTranspose.NTUmmfs <- function(mmfs)
{
  ret = list(n = mmfs$m,
             m = mmfs$n,
             nbX = mmfs$nbY,
             nbY = mmfs$nbX,
             neededNorm = normalizationTranspose(mmfs$neededNorm),
             A = t(mmfs$B),
             B = t(mmfs$A)
  )
  class(ret)="NTUmmfs"
  return(ret)
}
#
M.NTUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = mux0s[xs] * mmfs$A[xs,ys]
  term_2 = t( mu0ys[ys] * t(mmfs$B[xs,ys] ))
  #
  ret = pmin(term_1, term_2)
  #
  return(ret)
}
# 
dmux0s_M.NTUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mux0s * mmfs$A
  term_2 = t( mu0ys * t(mmfs$B ))
  #
  return(ifelse(term1 <= term2,1,0) * mmfs$A)
}
#
dmu0ys_M.NTUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mux0s * mmfs$A
  term_2 = t( mu0ys * t(mmfs$B ))
  #
  return(ifelse(term1 >= term2,1,0) * mmfs$B)  
}
#
dtheta_M.NTUmmfs <- function (mmfs,mux0s,mu0ys,dtheta=NULL)
{
  term_1 = mux0s * mmfs$A
  term_2 = t( mu0ys * t(mmfs$B ))
  t1lessthant2 = ifelse(term1 <= term2,1,0)
  der1 = mux0s * t1lessthant2
  der2 = t(mu0ys * t(1-t1lessthant2) )
  ret <- 0
  if(is.null(dtheta)){
    ret = cbind(Diagonal(x=der1),
                Diagonal(x=der2 ))
  }else{
    dtheta1 = matrix(dtheta[1:(mmfs$nbX*mmfs$nbY)],nrow = nbX)
    dtheta2 = matrix(dtheta[(1+mmfs$nbX*mmfs$nbY):(2*mmfs$nbX*mmfs$nbY)],nrow = nbX)
    ret = c(dtheta1 * der1 + dtheta2 * der2 )
  }
  return(ret)
}
#
margxInv.NTUmmfs = function(xs,mmfs,theMu0ys)
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  theaxs = inversePWA(mmfs$n[xs], t ( t(mmfs$B[xs,] / mmfs$A[xs,])* theMu0ys ),mmfs$A[xs,]) 
  return(theaxs)
}
#
margyInv.NTUmmfs = function(ys,mmfs,theMux0s)
{
  if (is.null(ys)) {ys = 1:length(mmfs$m)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  thebys = inversePWA(mmfs$m[ys], t(( mmfs$A[,ys] / mmfs$B[,ys]) * theMux0s), t(mmfs$B[,ys] ) ) 
  return(thebys)
}

################################################################################
########################            LTU MMfs            ########################
################################################################################
build_LTUmmfs <- function(n,m,lambda,K,neededNorm)
{
  ret = list(n=n,
             m=m,
             nbX=length(n),
             nbY=length(m),
             neededNorm=neededNorm,
             lambda = lambda,
             K=K,
             aux_zeta = 1-lambda)
  class(ret)="LTUmmfs"
  return(ret)
}
#
mmfsTranspose.LTUmmfs <- function(mmfs)
{
  ret = list(n=mmfs$m,
             m=mmfs$n,
             nbX = mmfs$nbY,
             nbY = mmfs$nbX,
             neededNorm=normalizationTranspose(mmfs$neededNorm),
             lambda = t(mmfs$aux_zeta),
             K=t(mmfs$K),
             aux_zeta = t(mmfs$lambda)
  )
  class(ret)="LTUmmfs"
  return(ret)
}
#
M.LTUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = mux0s[xs]^mmfs$lambda[xs,ys]
  term_2 = t( mu0ys[ys]^t(mmfs$aux_zeta[xs,ys]) )
  term_3 = mmfs$K[xs,ys]
  #
  ret = term_1 * term_2 * term_3
  #
  return(ret)
}
#
# 
dmux0s_M.LTUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mmfs$lambda* mux0s^(mmfs$lambda-1)
  term_2 = t( mu0ys^t(mmfs$aux_zeta) )
  term_3 = mmfs$K
  #
  ret = term_1 * term_2 * term_3
  #
  return(ret)
}
#
dmu0ys_M.LTUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mux0s^mmfs$lambda
  term_2 = mmfs$aux_zeta * t( mu0ys^t( mmfs$aux_zeta -1) )
  term_3 = mmfs$K
  #
  ret = term_1 * term_2 * term_3
  #
  return(ret)
}
#
dtheta_M.LTUmmfs <- function (mmfs,mux0s,mu0ys,dtheta=NULL)
{
  term_1 = mux0s^(mmfs$lambda)
  term_2 = t( mu0ys^t(mmfs$aux_zeta) )
  term_3 = mmfs$K
  logratio = log( mux0s %*% t(1 / mu0ys))
  der1 = logratio * term_1*term_2*term_3
  der2 = term_1 * term_2
  
  ret <- 0
  if(is.null(dtheta)){
    ret = cbind (Diagonal(x = der1),
                 Diagonal(x = der2))
  }else{
    dtheta1 = matrix(dtheta[1:(mmfs$nbX*mmfs$nbY)],nrow = nbX)
    dtheta2 = matrix(dtheta[(1+mmfs$nbX*mmfs$nbY):(2*mmfs$nbX*mmfs$nbY)],nrow = nbX)
    ret = c(dtheta1 * der1 + dtheta2 * der2)}
    
    return(ret)
}

#
################################################################################
########################            ETU MMfs            ########################
################################################################################
build_ETUmmfs <- function(n,m,C,D,kappa,neededNorm)
{
  ret = list(n=n,
             m=m,
             nbX=length(n),
             nbY=length(m),
             neededNorm=neededNorm,
             C = C,
             D = D,
             kappa = kappa
  )
  class(ret)="ETUmmfs"
  return(ret)
}
#
mmfsTranspose.ETUmmfs <- function(mmfs)
{
  ret = list(n = mmfs$m,
             m = mmfs$n,
             nbX = mmfs$nbY,
             nbY = mmfs$nbX,
             neededNorm = normalizationTranspose(mmfs$neededNorm),
             C = t(mmfs$D),
             D =  t(mmfs$C),
             tauinv = t(mmfs$tauinv)
  )
  class(ret)="ETUmmfs"
  return(ret)
}
#
M.ETUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = mmfs$C[xs,ys] * (mux0s[xs]^mmfs$kappa[xs,ys])
  term_2 = mmfs$D[xs,ys] * t(mu0ys[ys]^t(mmfs$kappa[xs,ys]))
  #
  ret = ((term_1 + term_2)/2)^(1/mmfs$kappa[xs,ys])
  #
  return(ret)
}
#
dmux0s_M.ETUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mmfs$C * (mux0s^mmfs$kappa)
  term_2 = mmfs$D * t(mu0ys^t(mmfs$kappa))
  #
  mus = ((term_1 + term_2)/2)^(1/mmfs$kappa)
  #
  
  return( mmfs$kappa * mmfs$C * (mus / mux0s)^(1 - mmfs$kappa) /2 )
}
#
dmu0ys_M.ETUmmfs <- function (mmfs,mux0s,mu0ys)
{
  term_1 = mmfs$C * (mux0s^mmfs$kappa)
  term_2 = mmfs$D * t(mu0ys^t(mmfs$kappa))
  #
  mus = ((term_1 + term_2)/2)^(1/mmfs$kappa)
  #
  
  return( mmfs$kappa * mmfs$D * t( (t(mus) / mu0ys)^t(1 - mmfs$kappa) ) /2 )
}
#

dtheta_M.ETUmmfs <- function (mmfs,mux0s,mu0ys,dtheta=NULL)
{
  term_1 = mmfs$C * (mux0s^mmfs$kappa)
  term_2 = mmfs$D * t(mu0ys^t(mmfs$kappa))
  #
  mus = ((term_1 + term_2)/2)^(1/mmfs$kappa)
  
  num = (mus / mmfs$kappa) * ( mmfs$C * mux0s^mmfs$kappa * log(mux0s) + mmfs$D * t(mu0ys^t(mmfs$kappa)*log(mu0ys)) - log(mus) )
  denom = mmfs$C * mux0s^mmfs$kappa  + mmfs$D * t(mu0ys^t(mmfs$kappa))
  
  der1 = mus^(1 - mmfs$kappa) * mux0s^mmfs$kappa /2
  der2 = mus^(1 - mmfs$kappa) * t(mu0ys^t(mmfs$kappa))/ 2
  der3 = num / denom
  
  ret <- 0
  if(is.null(dtheta)){
    ret = cbind (Diagonal(x = der1),
                 Diagonal(x = der2),
                 Diagonal(x = der3 ))
  }else{
    dtheta1 = matrix(dtheta[1:(mmfs$nbX*mmfs$nbY)],nrow = nbX)
    dtheta2 = matrix(dtheta[(1+mmfs$nbX*mmfs$nbY):(2*mmfs$nbX*mmfs$nbY)],nrow = nbX)
    dtheta3 = matrix(dtheta[(1+2*mmfs$nbX*mmfs$nbY):(3*mmfs$nbX*mmfs$nbY)],nrow = nbX)
    
    ret = c(dtheta1 * der1 + dtheta2 * der2 + dtheta3 * der3 )
  }
  
  return(ret)
}