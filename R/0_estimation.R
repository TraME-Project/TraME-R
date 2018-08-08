################################################################################
##
##   Copyright (C) 2015 - 2018 Alfred Galichon
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
dtheta_mu.DSE_model <- function(model, theta, deltatheta=diag(length(theta)))
{ 
  
  market = model$parametricMarket(model,theta)
  ### eventually, this will have to go
  ### ---------------- from here ------------------->>>>>>
  check_1 = (class(market$arumsG)=="logit") 
  check_2 = (class(market$arumsH)=="logit")
  check_3 = (market$arumsG$sigma == market$arumsH$sigma)
  if(check_1 && check_2 && check_3){
    ret = dtheta_mu_logit(model,market,theta,deltatheta)
    warning("Models of class DSE_model involve a market of type ITU-logit. This is inefficient and should be programmed using a market of type MFE.")
    return(ret)
  }
  ### <<<<<<---------- to here ---------------------------
  #
  outcome = solveEquilibrium(market,notifications=FALSE)
  #
  mu = outcome$mu
  mux0s = market$n-apply(mu,1,sum)
  mu0ys = market$m-apply(mu,2,sum)
  # 
  dtheta_Psi = model$dtheta_Psi
  dtheta_G = model$dtheta_G
  dtheta_H = model$dtheta_H
  
  U = outcome$U
  V = outcome$V
  

  deltaparamsPsi = dtheta_Psi(model, deltatheta=deltatheta)
  deltaparamsG = dtheta_G(model, deltatheta=deltatheta)
  deltaparamsH = dtheta_H(model, deltatheta=deltatheta)
  
  tr = market$transfers
  arumsG = market$arumsG
  arumsH = market$arumsH
  
  duPsimat = du_Psi(tr,U,V)
  dvPsimat = 1 - duPsimat
  duPsivec = c(duPsimat)
  dvPsivec = c(dvPsimat)
  
  #   duPsi = Diagonal(x = duPsimat) 
  #   dvPsi = Diagonal( x = dvPsimat)
  #
  HessGstar = D2Gstar(market$arumsG,mu,market$n,xFirst=T)  
  HessHstar = D2Gstar(market$arumsH,t(mu),market$m,xFirst=F)
  #
  denom = duPsivec * HessGstar + dvPsivec * HessHstar
  
  num1 = dparams_Psi(tr,U,V,deltaparamsPsi)
  num2 = duPsivec * dparams_NablaGstar(arumsG,mu,market$n,deltaparamsG,xFirst=TRUE)
  num3 = dvPsivec * dparams_NablaGstar(arumsH,t(mu),market$m,deltaparamsH,xFirst=FALSE)
  if (length(num1) == 0) (num1 = 0)
  if (length(num2) == 0) (num2 = 0)
  if (length(num3) == 0) (num3 = 0)
  #
  dmu = -solve(denom, num1+num2+num3)
  #num = cbind(num1, num2, num3)
  #
  return(list(mu= c(mu),mux0s=mux0s, mu0ys=mu0ys,dmu=dmu))
}
#
#
dtheta_mu.MFE_model <- function(model, theta, deltatheta=diag(length(theta)))
{
  
  market = model$parametricMarket(model,theta)
  if (class( market) != "MFE")
  {
    stop("MFE model does not involve a MFE market")
  }
  
  nbX = length(market$n)
  nbY = length(market$m)
  dtheta_params = model$dtheta_params
  #
  deltatheta = matrix(deltatheta , nrow = length(theta))
  rangeTheta = dim(deltatheta)[2]
  deltaparamsM = dtheta_params(model, theta=theta, deltatheta=deltatheta)
  
  
  mmfs = market$mmfs
  
  outcome = solveEquilibrium(market,notifications=FALSE,debugmode=FALSE)
  mu = outcome$mu
  mux0s = outcome$mux0
  mu0ys = outcome$mu0y
  
  du_Ms = matrix(dmux0s_M(mmfs,mux0s,mu0ys),nrow=nbX)
  dv_Ms = matrix(dmu0ys_M(mmfs,mux0s,mu0ys),nrow=nbX)
  # 
  deltaMs = matrix(dparams_M(mmfs,mux0s,mu0ys,deltaparamsM),nrow=nbX*nbY)
  
  deltaMs_array = array(deltaMs ,dim=c(nbX,nbY,rangeTheta))
  
  d_1 = apply(deltaMs_array, c(1,3), sum) / sigma
  d_2 = apply(deltaMs_array, c(2,3), sum) / sigma
  num =  - rbind(d_1,d_2)
  
  Delta11 = diag(1 + apply(du_Ms,1,sum),nrow=nbX)
  Delta22 = diag(1 + apply(dv_Ms,2,sum),nrow=nbY)
  Delta12 = dv_Ms
  Delta21 = t(du_Ms)
  Delta = rbind(cbind(Delta11,Delta12),cbind(Delta21,Delta22))
  
  dmusingles = solve(Delta,num)
  dmux0 = dmusingles[1:nbX,,drop=FALSE]
  dmu0y = dmusingles[(nbX+1):(nbX+nbY),,drop=FALSE]
  dmux0full = array(0,dim=c(nbX,nbY,rangeTheta))
  dmu0yfull = array(0,dim=c(nbX,nbY,rangeTheta))
  
  for(y in 1:nbY){
    dmux0full[,y,] = dmux0
  }
  for(x in 1:nbX){
    dmu0yfull[x,,] = dmu0y
  }
  
  dmu = c(du_Ms)*matrix(dmux0full, ncol=rangeTheta) + 
    c(dv_Ms)*matrix(dmu0yfull, ncol=rangeTheta) +
    deltaMs     
  
  return(list(mu = c(mu), mux0s = mux0s, mu0ys = mu0ys, dmu = dmu))
  
}
#
#
#

mLogLikelihood <- function(theta, model, muhat, muhatx0, muhat0y, scale=1, byIndiv=T) # to be modified
{
  mudmu = try( dtheta_mu(model,theta),silent=T)
  #
  ret <- 0
  if(class(mudmu)!="try-error"){
    if (byIndiv==TRUE){
      mLL = - sum(2* muhat * log(mudmu$mu)) - sum(muhatx0 * log(mudmu$mux0s)) - sum(muhat0y * log(mudmu$mu0ys))
      #
      term_1 = t(2*muhat/matrix(mudmu$mu,nrow=model$nbX) - muhatx0/mudmu$mux0s)
      term_2 = muhat0y / mudmu$mu0ys
      term_grad = c(t(term_1 - term_2))*mudmu$dmu
      #
      mGradLL = - apply(term_grad,2,sum) 
      
      
    } else {
      N = sum(c(c(mudmu$mu),c(mudmu$mux0s), c(mudmu$mu0ys)))
      mLL = - sum(muhat * log(mudmu$mu/N)) - sum(muhatx0 * log(mudmu$mux0s/N)) - sum(muhat0y * log(mudmu$mu0ys/N))
      #
      term_1 = t(muhat/matrix(mudmu$mu,nrow=model$nbX) - muhatx0/mudmu$mux0s)
      term_2 = muhat0y / mudmu$mu0ys
      term_grad = c(t(term_1 - term_2))*mudmu$dmu
      #
      mGradLL = - apply(term_grad,2,sum)
      term_3 = sum(c(muhat, muhatx0, muhat0y))/N * apply(mudmu$dmu,2,sum)
      mGradLL = mGradLL - term_3
      
    }
    #
    ret = list(objective = mLL / scale, gradient = mGradLL / scale)
  }else{
    ret = list(objective=NaN, gradient=rep(NaN,model$dimTheta))
  }
  #
  return(ret)
}

mle <- function(model, muhat, theta0=NULL, xtol_rel=1e-8, maxeval=1e5, print_level=0, byIndiv=T)
{
  inittheta = model$inittheta
  nbX = length(model$n)
  nbY = length(model$m)
  scale = max(sum(model$n),sum(model$n))
  #dimTheta = length(model$dimTheta) #Keith: should this be model$dimTheta?
  dimTheta = model$dimTheta
  if(print_level > 0){
    message(paste0("Maximum Likelihood Estimation of ",class(model)," model."))
  }
  #
  muhatx0  = model$n-apply(muhat,1,sum)
  muhat0y  = model$m-apply(muhat,2,sum)
  #
  if(is.null(theta0)){
    theta0 = inittheta(model)$theta
  }
  #
  lb     = inittheta(model)$lb
  ub     = inittheta(model)$ub
  #
  tm = proc.time()
  res = nloptr(x0=theta0, 
               eval_f=mLogLikelihood,
               lb=lb, ub=ub,
               opt = list(algorithm='NLOPT_LD_LBFGS',
                          xtol_rel=xtol_rel,
                          maxeval=maxeval,
                          "print_level"=print_level), 
               model=model,
               muhat=muhat,
               muhatx0=muhatx0,
               muhat0y=muhat0y,
               scale = scale,
               byIndiv=byIndiv)
  time = proc.time() - tm
  time = time["elapsed"]
  #
  if(print_level > 0){
    print(res, show.controls=((1+nbX*nbY):(dimTheta+nbX*nbY)))
  }
  #
  return(list(thetahat=res$solution, fval = res$objective, time=time))
}


plotCovariogram2D = function(model,lambda,muhat = NULL, dim1=1,dim2=2,nbSteps =100)
{
  # HERE, INSERT FILTER TO APPLY ONLY TO TU MODELS W LINEAR PARAMETERIZATION
  library(gplots)
  points = matrix(0,nbSteps,2)
  for (i in (1:nbSteps))
  {
    angle=2 * pi * i / (nbSteps-1)
    lambdastar = lambda
    lambdastar[dim1]=lambda[dim1]*cos(angle) - lambda[dim2]*sin(angle) 
    lambdastar[dim2]=lambda[dim1]*sin(angle) + lambda[dim2]*cos(angle) 
    market = model$parametricMarket(model,lambdastar)
    mustar = solveEquilibrium(market)$mu
    points[i,] =  c(c(mustar) %*% matrix(phi_xyk, ncol=model$dimTheta) ) [c(dim1,dim2)]
    
  }
  par(mar = rep(2, 4))
  plot(points,col=rgb(0, 0,1))
  lines(points,col=rgb(0, 0,1))
  
  if (! is.null(muhat))
  {
    Chat = c(c(muhat) %*% matrix(phi_xyk, ncol=model$dimTheta) ) [c(dim1,dim2)]
    points(matrix(Chat,1,2),pch=16,col=rgb(1,0,0))
  }
}
