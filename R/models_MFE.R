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
################################################################################
########################       default methods           #######################
################################################################################
inittheta_default <- function(model)
{
  ret = list(theta=rep(0,model$dimTheta),
             lb = NULL, ub = NULL)
}
#
################################################################################
########################       affinity model            #######################
################################################################################
#
#
buildModel_affinity <- function(Xvals, Yvals, n=NULL, m=NULL, sigma = 1 )
{
  nbX = dim(Xvals)[1]
  nbY = dim(Yvals)[1]
  #
  dX = dim(Xvals)[2]
  dY = dim(Yvals)[2]
  #
  if(is.null(n)){
    n = rep(1,nbX)
  }
  if(is.null(m)){
    m = rep(1,nbY)
  }
  #
  if ( sum(n) != sum(m) ) {stop("Unequal mass of individuals in an affinity model.")}
  #
  neededNorm = defaultNorm(TRUE)
  #
  ret = list(type = c("affinity"),
             dimTheta=dX*dY,
             dX=dX, dY=dY,
             nbX = nbX, nbY = nbY,
             n=n, m=m,
             sigma = sigma,
             neededNorm = neededNorm,
             phi_xyk_aux = kronecker(Yvals,Xvals),
             Phi_xyk = function(model)
               (model$phi_xyk_aux),
             Phi_xy = function(model ,lambda)
               ( array(model$phi_xyk_aux,dim=c(model$nbX*model$nbY,model$dimTheta)) %*% lambda ),
             Phi_k = function(model, muhat)
               (c(c(muhat) %*% array(model$phi_xyk_aux,dim=c(model$nbX*model$nbY,model$dimTheta)))),
             inittheta = inittheta_default,
            # dtheta_Psi = dtheta_Psi_affinity,
            dtheta_params = dtheta_params_affinity,
            # dtheta_G = dtheta_G_default,
            # dtheta_H = dtheta_H_default,
             mme = mme_affinity,
             parametricMarket = parametricMarket_affinity
  )

  class(ret) = "MFE_model"
  #
  return(ret)
}
#
dtheta_params_affinity <- function(model, theta, deltatheta=diag(model$dimTheta))
  (return( exp(  matrix(model$Phi_xy(model,c(theta)), nrow=model$nbX)/2) * model$Phi_xy(model, deltatheta) /2 ))
#
#
mmeaffinityNoRegul <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
  # mmeaffinityNoRegul should be improved as one should make use of the logit structure and use the ipfp
{
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of Affinity model via IPFP+BFGS optimization."))
  }
  #
  theta0 = model$inittheta(model)$theta
  Chat = model$Phi_k(model, muhat)
  nbX = model$nbX
  nbY = model$nbY
  dX = model$dX
  dY = model$dY
  dimTheta = model$dimTheta
  sigma = model$sigma
  #
  totmass = sum(model$n)

  if ( sum(model$n) != totmass ) {stop("Unequal mass of individuals in an affinity model.")}
  if (sum(muhat) !=  totmass) {stop("Total number of couples does not coincide with margins.")}
  p = model$n / totmass
  q = model$m / totmass
  IX=rep(1,nbX)
  tIY=matrix(rep(1,nbY),nrow=1)
  f = p %*% tIY
  g = IX %*% t(q)
  pihat = muhat / totmass
  v=rep(0,nbY)
  #
  #
  valuef = function(A)
  {
    Phi = matrix(model$Phi_xy(model,c(A)) ,nbX,nbY)
    # Phi = Xvals %*% matrix(A,nrow=dX) %*% t(Yvals)
    contIpfp = TRUE
    iterIpfp = 0
    while(contIpfp)
    {
      iterIpfp = iterIpfp+1
      u = sigma*log(apply(g * exp( ( Phi - IX %*% t(v) ) / sigma ),1,sum))
      vnext = sigma*log(apply(f * exp( ( Phi - u %*% tIY ) / sigma ),2,sum))
      error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
      if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
      v=vnext
    }
    #print(c("Converged in ", iterIpfp, " iterations."))
    pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / model$sigma )
    if (iterIpfp >= maxiterIpfp ) {stop('maximum number of iterations reached')}
    v <<- vnext
    #thegrad =  c(    c(pi - pihat) %*% phis)
    thegrad = model$Phi_k(model,pi - pihat)
    theval = sum(thegrad * c(A)) - sigma* sum(pi*log(pi))

    return(list(objective = theval,gradient = thegrad))
  }

  A0 = rep(0,dX*dY)
  resopt = nloptr(x0=A0,
                  eval_f = valuef,
                  opt = list(algorithm = 'NLOPT_LD_LBFGS',
                             xtol_rel = xtol_rel,
                             maxeval=maxeval,
                             print_level = print_level))
  #  AffinityMatrix = matrix(res$solution,nrow=dX)
  if (resopt$status<0) {warning("nloptr convergence failed.")}
  #
  thetahat = resopt$solution
  ret =list(thetahat=thetahat,
            val=resopt$objective,
            status = resopt$status)
  #
  return(ret)
}
#
#
mmeaffinityWithRegul <- function(model, muhat, lambda, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
  # Reference: Arnaud Dupuy, Alfred Galichon, Yifei Sun (2016). "Learning Optimal Transport Costs under Low-Rank Constraints."
  # Implementation by Yifei Sun.
{
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of Affinity model with regularization via proximal gradient."))
  }
  #
  theta0 = model$inittheta(model)$theta
  Chat = model$Phi_k(model, muhat)
  nbX = model$nbX
  nbY = model$nbY
  dX = model$dX
  dY = model$dY
  dimTheta = model$dimTheta
  sigma = model$sigma
  #
  totmass = sum(model$n)
  #
  if ( sum(model$n) != totmass ) {stop("Unequal mass of individuals in an affinity model.")}
  if (sum(muhat) !=  totmass) {stop("Total number of couples does not conicide with margins.")}
  p = model$n / totmass
  q = model$m / totmass
  IX=rep(1,nbX)
  tIY=matrix(rep(1,nbY),nrow=1)
  f = p %*% tIY
  g = IX %*% t(q)
  pihat = muhat / totmass
  v=rep(0,nbY)
  A = rep(0,dX*dY)
  t_k = .3   # step size for the prox grad algorithm (or grad descent when lambda=0)
  iterCount = 0
  while (1)
  {
    # compute pi^A
    # Phi = Xvals %*% matrix(A,nrow=dX) %*% t(Yvals)
    Phi = matrix(model$Phi_xy(model,c(A)) ,nbX,nbY)
    contIpfp = TRUE
    iterIpfp = 0
    while(contIpfp)
    {
      iterIpfp = iterIpfp+1
      u = sigma*log(apply(g * exp( ( Phi - IX %*% t(v) ) / sigma ),1,sum))
      vnext = sigma*log(apply(f * exp( ( Phi - u %*% tIY ) / sigma ),2,sum))
      error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
      if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
      v=vnext
    }

    pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / sigma )

    if (iterIpfp >= maxiterIpfp ) {stop('maximum number of iterations reached')}

    # do prox grad descent
    # thegrad = c(phis %*% c(pi - pihat))
    thegrad = model$Phi_k(model,pi - pihat)


    # take one gradient step
    A = A - t_k*thegrad

    if (lambda > 0)
    {
      # compute the proximal operator
      SVD = svd(matrix(A,nrow=dX))
      U = SVD$u
      D = SVD$d
      V = SVD$v

      D = pmax(D - lambda*t_k, 0.0)
      A = c(U %*% diag(D) %*% t(V))
    } # if lambda = 0 then we are just taking one step of gradient descent
    ### testing optimality
    if (iterCount %% 10 == 0)
    {
      alpha = 1.0
      tmp = svd(matrix(A - alpha * thegrad, nrow=dX))
      tmp_second = sum((A - c(tmp$u %*% diag(pmax(tmp$d - alpha*lambda, 0.0)) %*% t(tmp$v)))^2)
      cat("testing optimality ", tmp_second, "\n")
    }

    if (lambda > 0)
    {
      theval = sum(thegrad * c(A)) - sigma * sum(pi*log(pi)) + lambda * sum(D)
    } else
    {
      theval = sum(thegrad * c(A)) - sigma * sum(pi*log(pi))
    }

    iterCount = iterCount + 1

    if (iterCount>1 && abs(theval - theval_old) < 1E-6) { break }

    theval_old = theval

  }
  ret =list(thetahat=c(A),
            val=theval)

  #
  return(ret)
}
#
#
mme_affinity <- function(model, muhat, lambda = NULL, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
{
  if (is.null(lambda)) {lambda = 0}
  if ( lambda == 0 )
  {return(mmeaffinityNoRegul(model,muhat, xtol_rel , maxeval , tolIpfp , maxiterIpfp , print_level))}
  else
  { return(mmeaffinityWithRegul(model,muhat, lambda, xtol_rel , maxeval , tolIpfp , maxiterIpfp , print_level)) }
}
#
parametricMarket_affinity <- function(model, theta)
  (build_market_geoMFE(model$n,model$m,
                       exp(  matrix(model$Phi_xy(model,c(theta)), nrow=model$nbX)/2),
                       neededNorm=model$neededNorm))
