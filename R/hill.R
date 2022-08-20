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



#' @title Estimation of the Half Maximal Inhibitory Concentration
#' 
#' @description The half maximal inhibitory concentration (IC50) is a measure of the
#' effectiveness of a compound in inhibiting biological or biochemical
#' function.  This quantitative measure indicates how much of a particular drug
#' or other substance (inhibitor) is needed to inhibit a given biological
#' process (or component of a process) by half.
#' 
#' See reference(s).
#' 
#' In this version of the function the maximal occupancy (rmax) is estimated
#' automatically.  This should be optional.
#' 
#' @param conc a vector of drug concentrations in plasma (example units are
#' ng/mL).
#' @param occ a vector of PET occupancy values that correspond to the measured
#' drug concentrations in plasma.
#' @param guess a length-two vector of starting values for the nonlinear
#' optimization.
#' @param control is a list of parameters used by \code{nls.lm.control} that
#' are set by default, but may be customized by the user.
#' @return List with the following elements \itemize{
#' \item{IC50}{Half maximal inhibitory concentration}
#' \item{rmax}{Estimated maximal occupancy} \item{IC50SE}{Approximate standard
#' error for IC50} \item{rmaxSE}{Approximate standard erorr for rmax}
#' \item{hessian}{Hessian matrix from the Levenburg-Marquardt procedure}
#' \item{info}{Return value from the Levenburg-Marquardt procedure}
#' \item{deviance}{Deviance from the Levenburg-Marquardt procedure}
#' \item{message}{Text message from the Levenburg-Marquardt procedure}
#' }
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link[minpack.lm]{nls.lm}}
#' @references \href{https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)}{Hill
#' Equation}
#' \href{https://en.wikipedia.org/wiki/IC50}{IC50}
#' @export hillEquation
#' @importFrom stats coef median
hillEquation <- function(conc, occ, guess=c(1,100), control=minpack.lm::nls.lm.control()) {
    func <- function(x, conc, occ) {
    IC50 <- x[1]
    rmax <- x[2]
    occ - rmax * conc / (conc + IC50)
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control, conc=conc, occ=occ)
  rdf <- length(out$fvec) - length(coef(out))
  varcovmat <- (out$deviance / rdf) * chol2inv(chol(out$hessian))
  list(IC50=out$par[1],
       rmax=out$par[2],
       IC50SE=sqrt(varcovmat[1,1]),
       rmaxSE=sqrt(varcovmat[2,2]),
       hessian=out$hessian, info=out$info, deviance=out$deviance,
       message=out$message)
}
