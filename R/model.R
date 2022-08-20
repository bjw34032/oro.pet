##
##
## Copyright (c) 2010 Brandon Whitcher and Volker Schmid
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
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
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
##
## $Id:$
##



#' Compartmental Models for Kinetic Parameter Estimation
#' 
#' A selection of parametric models are provided that combine a compartmental
#' model for tissue and empirical versions of the arterial input function or
#' reference region time activity curve.
#' 
#' Parametric models from the PET literature are provided to the user for
#' kinetic parameter estimation.
#' 
#' @param type is a character string that identifies the type of compartmental
#' model to be used.  Acceptable models include: \describe{
#' \item{list("srtm")}{Simplified Reference Tissue Model}
#' \item{list("srtm2")}{Simplified Reference Tissue Model in two steps} }
#' @return A function.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{simplifiedReferenceTissueModel}}
#' @references
#' 
#' Lammertsma, A.A and Hume, S.P. (1996) Simplified reference tissue model for
#' PET receptor studies, \emph{NeuroImage}, \bold{4}, 153-158.
#' 
#' Wu, Y and Carson, R.E. (2002) Noise reduction in the simplified reference
#' tissue model for neuroreceptor functional imaging, \emph{Journal of Cerebral
#' Blood Flow & Metabolism}, \bold{22}, 1440-1452.
#' @keywords misc
#' @export compartmentalModel
#' @importFrom stats approx
compartmentalModel <- function(type) {
  switch(type,
         srtm =
         function(time, theta, ref, parameter.not.used) {
           R1 <- theta[1]
           k2prime <- theta[2]
           k2 <- theta[3]
           tsec <- seq(min(time * 60), ceiling(max(time * 60)), by=1)
           tsec.gt0 <- tsec[tsec >= 0]
           ltsec.gt0 <- length(tsec.gt0)
           tsec.lt0 <- tsec[tsec < 0]
           ref.sec <- approx(time * 60, ref, tsec.gt0)$y
           if (is.na(ref.sec[ltsec.gt0])) {
             ref.sec[ltsec.gt0] <- ref.sec[ltsec.gt0 - 1]
           }
           erg.gt0 <- approx(tsec.gt0,
                             expConv(ref.sec,
                                     R1 * (k2prime - k2) / 60,
                                     k2 / 60),
                             time[time >= 0] * 60)$y
           erg.gt0[1] <- 0
           erg <- numeric(length(time))
           erg[time >= 0] <- R1 * ref[time >= 0] + erg.gt0
           return(erg)
         },
         srtm2 =
         function(time, theta, ref, k2prime) {
           th1 <- theta[1]
           th2 <- theta[1] * (k2prime - theta[2])
           th3 <- theta[2]
           tsec <- seq(min(time * 60), ceiling(max(time * 60)), by=1)
           tsec.gt0 <- tsec[tsec >= 0]
           ltsec.gt0 <- length(tsec.gt0)
           tsec.lt0 <- tsec[tsec < 0]
           ref.sec <- approx(time * 60, ref, tsec.gt0)$y
           if (is.na(ref.sec[ltsec.gt0])) {
             ref.sec[ltsec.gt0] <- ref.sec[ltsec.gt0 - 1]
           }
           erg.gt0 <- approx(tsec.gt0,
                             expConv(ref.sec, th2 / 60, th3 / 60),
                             time[time >= 0] * 60)$y
           erg.gt0[1] <- 0
           erg <- numeric(length(time))
           erg[time >= 0] <- th1 * ref[time >= 0] + erg.gt0
           return(erg)
         })
}



#' Empirical Convolution Between an Input Function and a Single Exponential
#' 
#' Computationally efficient method to convolve a vector of observations and a
#' single exponential function with two parameters.
#' 
#' Assuming the input function has been sampled (or interpolated) to a high
#' temporal resolutions, say one Hertz, a simple for loop is used to perform
#' the convolution.
#' 
#' @param input is the so-called input function.
#' @param k1 is the scaling parameter in the single exponential function.
#' @param k2 is the decay parameters in the single exponential function.
#' @return The vector containing the result from the convolution operation.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @keywords misc
#' @export expConv
expConv <- function(input, k1, k2) {
  k1input <- k1 * input
  if (k2 == 0) {
    convolution <- cumsum(k1input)
  } else {
    prev <- 0
    len <- length(input)
    convolution <- numeric(len)
    ek2 <- exp(-k2)
    k1intputk2 <- k1input * (1 - ek2) / k2
    for (i in 1:len) {
      prev <- prev * ek2 + k1intputk2[i]
      convolution[i] <- prev
    }
  }
  return(convolution)
}
  
