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
# GEO <-> TU
# CES <-> ETU
# CD <-> LTU
# MIN <-> NTU


