################################################################################
########################       affinity model            #######################
################################################################################
#
#
buildModel_affinityDSE <- function(Xvals, Yvals, n=NULL, m=NULL )
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
             dtheta_Psi = dtheta_Psi_affinity,
             dtheta_G = dtheta_G_default,
             dtheta_H = dtheta_H_default,
             mme = mme_affinity,
             parametricMarket = parametricMarket_affinity_DSE
  )
  
  class(ret) = "DSE_model"
  #
  return(ret)
}
#
dtheta_Psi_affinity <- function(model, theta=NULL, deltatheta=diag(model$dimTheta))
  (return(model$Phi_xy(model, deltatheta)))
#
#
parametricMarket_affinity_DSE <- function(model, theta) 
  (build_market_TU_logit(model$n,model$m,
                         matrix(model$Phi_xy(model,c(theta)), nrow=model$nbX),
                         sigma=model$sigma,neededNorm=model$neededNorm))

