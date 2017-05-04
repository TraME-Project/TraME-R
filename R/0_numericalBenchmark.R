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


dtheta_mu_numeric <- function (model, theta, deltatheta=diag(length(theta)))
{ 
  market = model$parametricMarket(model,theta)
  thef <- function(thetheta){
    ret = solveEquilibrium(model$parametricMarket(model,thetheta),notifications=FALSE)$mu
    return(ret)
  }
  #
  outcome = solveEquilibrium( model$parametricMarket(model,theta),notifications=FALSE)
  #
  mu=outcome$mu
  mux0s = outcome$mux0
  mu0ys = outcome$mu0y
  #
  dmu = jacobian(thef,theta) %*% deltatheta
  #
  return(list(mu = c(mu), mux0s = mux0s, mu0ys = mu0ys, dmu = dmu))
}