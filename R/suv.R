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

#' @rdname suv
#' @title Calculating the Lean Body Mass
#' 
#' @description The lean body mass (LBM) is calculated according to the formula
#' \deqn{1.1\cdot\mbox{weight}-128\cdot(\mbox{weight}/\mbox{height})^2} if male
#' and \deqn{1.07\cdot\mbox{weight}-148\cdot(\mbox{weight}/\mbox{height})^2} if
#' female.
#' 
#' 
#' @aliases leanBodyMass
#' @param height is a vector of heights in centimeters.
#' @param weight is a vector of weights in kilograms.
#' @param gender is a character vector (may be of length one) with the value
#' \dQuote{male} or \dQuote{female}.
#' @return Vector of lean body mass values in kilograms.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{standardUptakeValue}}
#' @references Sugawara, Y., K. R. Zasadny, A. W. Neuhoff, R. L. Wahl (1999)
#' Reevaluation of the Standardized Uptake Value for FDG: Variations with Body
#' Weight and Methods for Correction, \emph{Radiology} \bold{213}: 521--525.
#' @examples
#' library(oro.pet)
#' n <- 11
#' h <- seq(200, 150, length=n)
#' w <- seq(80, 120, length=n)
#' cbind(h, w, leanBodyMass(h, w, "male"), leanBodyMass(h, w, "female"))
#' @export
leanBodyMass <- function(height, weight, gender) {
  ## weight in kg
  ## height in cm
  if (length(height) != length(weight)) {
    stop("Length of height and weight vectors must be equal")
  }
  n <- length(height)
  lbm <- numeric(n)
  if (length(gender) == 1) {
    gender <- rep(gender, n)
  } else {
    if (length(gender) != length(height)) {
      stop("Length of gender vector must equal height and weight")
    }
  }
  m <- gender == "male"
  lbm[m] <- 1.1 * weight[m] - 128 * (weight[m] / height[m])^2
  f <- gender == "female"
  lbm[f] <- 1.07 * weight[f] - 148 * (weight[f] / height[f])^2
  return(lbm)
}

#' @name Summarizing SUVs
#' @rdname suv
#' @title Summarizing SUVs for PET
#' 
#' The standard uptake value (SUV) is summarized using the hotspot method or by
#' calculating total volume of the high values.
#' 
#' 
#' @aliases hotSpotSUV totalSUV
#' @param suv is the standard uptake value (SUV).
#' @param radius is the desired hotspot radius (units = voxels).
#' @param type is a character string (acceptable values are \code{2D} or
#' \code{3D}) that determines the dimension of the hot spot (default =
#' \code{3D}).
#' @param mask is a multidimensional array of logical values.
#' @param z is the slice index.
#' @param bg is the estimated background SUV.
#' @param local is a logical value.
#' @return ...
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @export
#' @seealso \code{\link{leanBodyMass}}
#' @importFrom oro.nifti pixdim
#' @importFrom stats median
hotSpotSUV <- function(suv, radius=10, type="3D") {
  circle <- function(XY, center, r, pixdim=pixdim) {
    x <- matrix(1:XY[1], XY, byrow=TRUE) * pixdim[1]
    y <- matrix(1:XY[2], XY) * pixdim[2]
    center <- center * pixdim
    (center[1] - x)^2 + (center[2] - y)^2 <= r^2
  }
  sphere <- function(XYZ, center, r, pixdim=pixdim) {
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
    hotSpotMask <- circle(dim(suv)[1:2], m[2:1], radius, pixdim(suv)[2:3]) # Why m[2:1]?
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

#' @rdname suv
#' @export
totalSUV <- function(suv, mask, z, bg, local=TRUE) {
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

