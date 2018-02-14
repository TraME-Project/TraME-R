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
marketTranspose.MFE <- function(market)
{
  thelist=list(n=market$m, m=market$n, 
               neededNorm=normalizationTranspose(market$neededNorm),
               mmfs = mmfsTranspose(market$mmfs)
               )
  #
  add = 0
  names = names(market)
  # here, transpose additional elements; add 1 to add each time
  if(length(names) > length(thelist) + add){
    message("Warning: in bipartite market transposition (MFE), 
            some elements have not been copied.")
  }
  #
  return(structure(thelist,class=class(market)))
  }


#################################################
########     Methods for MFE markets    #########
#################################################
#
DSEToMFE <- function(market)
{
  if ((class(market) != "DSE")  | (class(market$arumsG) != "logit") | (class(market$arumsH) != "logit") | (market$arumsG$sigma != 1 ) | (market$arumsG$sigma != 1) )
  { stop("market object to convert is not DSE; or not logit with sigma=1") }
  mmfs = PsiToM(market$transfers,market$n,market$m,market$neededNorm)
  
  ret = list(kind = class(mmfs),
             n=market$n, m=market$m,
             neededNorm=market$neededNorm,
             mmfs=mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)  
}

# GEO <-> TU
build_market_geoMFE <- function(n, m, K, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  mmfs = build_geommfs(n,m,K,neededNorm)
  #
  ret = list(kind = class(mmfs),
             n=n,m=m,
             neededNorm=neededNorm,
             mmfs =mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)
}
# CES <-> ETU
build_market_cesMFE <- function(n, m, C,D,kappa, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  mmfs = build_cesmmfs(n,m,C,D,kappa,neededNorm)
  #
  ret = list(kind = class(mmfs),
             n=n,m=m,
             neededNorm=neededNorm,
             mmfs =mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)
}


# cod <-> LTU
build_market_codMFE <- function(n, m, lambda, K , neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  if(is.null(neededNorm)){
    outsideOption = TRUE
  }else{
    outsideOption = FALSE
  }
  mmfs = build_codmmfs(n,m,lambda,K,neededNorm)
  #
  ret = list(kind = class(mmfs),
             n=n,m=m,
             neededNorm=neededNorm,
             mmfs =mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)
}

# MIN <-> NTU
build_market_minMFE <- function(n, m, A, B, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  mmfs = build_minmmfs(n,m,A,B,neededNorm)
  #
  ret = list(kind = class(mmfs),
             n=n,m=m,
             neededNorm=neededNorm,
             mmfs =mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)
}

# general MFE markets
build_market_MFE <- function(mmfs)
{
  neededNorm=mmfs$neededNorm
  n=mmfs$n 
  m=mmfs$m 
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  if(is.null(neededNorm)){
    outsideOption = TRUE
  }else{
    outsideOption = FALSE
  }
  #
  ret = list(kind = class(mmfs),
             n=n,m=m,
             neededNorm=neededNorm,
             mmfs =mmfs
  )
  class(ret) = "MFE"
  #
  return(ret)
}
