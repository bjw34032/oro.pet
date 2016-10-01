##
## Copyright (c) 2012, Brandon Whitcher
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



#' Compute Drug Occupancy with Approximate Standard Errors
#' 
#' Receptor occupancy is calculated from posititron emission tomography (PET)
#' data as the treatment-induced relative change in the concentration of
#' available (not occupied) receptors.
#' 
#' Occupancy is calculated using the straightforward and well-known formula.
#' If the standard errors for the two binding potentials are provided, then the
#' delta method is used to approximate the standard error for the estimate of
#' occupancy.
#' 
#' @param base is the baseline binding potential (BPND).
#' @param drug is the post-treatment binding potential (BPND).
#' @param baseSE is the standard error for the baseline BPND.
#' @param drugSE is the standard error for the post-treatment BPND.
#' @param base.drug.corr is the user-specified correlation between baseline and
#' post-treatment binding potentials.
#' @return \item{OCC}{is the percent drug occupancy.} \item{SE}{is the
#' approximate standard error of the parameter estimate.}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link[msm]{deltamethod}}
#' @references Cunningham VJ, Rabiner EA, Slifstein M, Laruelle M (2010).
#' Measuring drug occupancy in the absence of a reference region: the Lassen
#' plot re-visited, \emph{Journal of Cerebral Blood Flow & Metababolism},
#' \bold{30}, 46-50.
#' 
#' Passchier J, Gee A, Willemsen A, Vaalburg W, van Waarde A (2002).  Measuring
#' drug-related receptor occupancy with positron emission tomography,
#' \emph{Methods}, \bold{27}, 278-286.
#' @export occupancy
occupancy <- function(base, drug, baseSE=NULL, drugSE=NULL,
                      base.drug.corr=0) {
  if (length(base) != length(drug)) {
    stop("Length of binding potential vectors must be equal")
  }
  occ.se <- rep(NA, length(base))
  if (! is.null(baseSE) && ! is.null(drugSE)) {
    if ((is.null(baseSE) && ! is.null(drugSE)) ||
        (! is.null(baseSE) && is.null(drugSE))) {
      stop("Both sets of standard errors must be provided")
    }
    if (length(baseSE) != length(drugSE)) {
      stop("Length of SE vectors must be equal")
    }
    for (k in 1:length(baseSE)) {
      x <- c(base[k], drug[k])
      x.cov <- base.drug.corr * baseSE[k] * drugSE[k]
      varcov <- matrix(c(baseSE[k]^2, x.cov, x.cov, drugSE[k]^2), 2, 2)
      occ.se[k] <- msm::deltamethod(~ (x1 - x2) / x1, x, varcov)
    }
  }
  list(OCC = (base - drug) / base, SE = occ.se)
}
