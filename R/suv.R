##
## Copyright (c) 2009-2011, Brandon Whitcher
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

leanBodyMass <- function(height, weight, gender="male") {
  ## weight in kg
  ## height in cm
  switch(gender,
         "male" =  1.1 * weight - 128 * (weight / height)^2,
         "female" = 1.07 * weight - 148 * (weight / height)^2,
         stop("Incorrect gender specified"))
}

standardUptakeValue <- function(data, mask, dose, mass) {
  ## Make sure the units cancel!
  out <- array(NA, dim(data))
  out[mask] <- data[mask] / dose * mass
  as(out, "nifti") <- data
  return(out)
}

hotSpotSUV <- function(suv, radius=10, type="2D") {
  circle <- function(XY, center, r) {
    x <- matrix(1:XY[1], XY, byrow=TRUE)
    y <- matrix(1:XY[2], XY)
    (center[1] - x)^2 + (center[2] - y)^2 <= r^2
  }
  sphere <- function(XYZ, center, r, pixdim=c(1,1,1)) {
    x <- array(1:XYZ[1], XYZ) * pixdim[1]
    y <- array(1:XYZ[2], XYZ[c(2,1,3)])
    y <- aperm(y, c(2,1,3)) * pixdim[2]
    z <- array(1:XYZ[3], XYZ[c(3,1,2)])
    z <- aperm(z, c(2,3,1)) * pixdim[3]
    center <- center * pixdim
    (center[1] - x)^2 + (center[2] - y)^2 + (center[3] - z)^2 <= r^2
  }
  suv.max <- max(suv, na.rm=TRUE)
  m <- which(suv == suv.max, arr.ind=TRUE)
  if (type == "2D") {
    hotSpotMask <- circle(dim(suv)[1:2], m[2:1], radius) # Why m[2:1]?
  } else {
    hotSpotMask <- sphere(dim(suv), m, radius, pixdim(suv)[2:4])
  }
  hotSpotMask <- ifelse(hotSpotMask > 0, TRUE, FALSE)
  list(mean = mean(suv[hotSpotMask], na.rm=TRUE),
       median = median(suv[hotSpotMask], na.rm=TRUE),
       min = min(suv[hotSpotMask], na.rm=TRUE),
       max = suv.max,
       voxels = suv[hotSpotMask],
       mask = hotSpotMask)
}

totalSUV <- function(suv, mask, z, bg, local=TRUE, mname=NULL, nt=NULL) {
  j <- which(suv == max(suv, na.rm=TRUE), arr.ind=TRUE)
  suv.max <- suv[j]
  Z <- length(z)
  vol <- numeric(Z)
  vol.mask <- array(FALSE, dim(suv))
  total <- vector("list", Z)
  for (i in 1:Z) {
    sli <- suv[, , z[i]]
    if (local) {
      k <- which(sli == max(sli, na.rm=TRUE), arr.ind=TRUE)
      sli.max <- sli[k]
      thresh <- mean(c(bg, sli.max))
    } else {
      thresh <- mean(c(bg, suv.max))
    }
    vol.mask[, , z[i]] <- sli >= thresh & mask[, , z[i]]
    vol.mask[is.na(vol.mask)] <- FALSE
    vol[i] <- sum(vol.mask[, , z[i]], na.rm=TRUE)
    total[[i]] <- sli[vol.mask[, , z[i]]]
  }
  list(volume = sum(vol), mean = mean(unlist(total)),
       median = median(unlist(total)), mask = vol.mask)
}

