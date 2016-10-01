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



#' Plot Baseline Versus Post-Treatment Binding Potentials
#' 
#' Inspired by the Lassen plot (Cunningham et al., 2010) this is a
#' straightforward graphical summary of pre-treatment versus post-treatment
#' binding potentials for a single subject across multiple brain regions.
#' 
#' See the reference below.
#' 
#' @param base is the vector of baseline binding potentials across brain
#' regions.
#' @param drug is the vector of post-treatment binding potentials across brain
#' regions.
#' @param lty45 is the line type for the 45-degree line.
#' @param lty is the line type for the estimated regression line.
#' @param lwd45 is the line width for the 45-degree line.
#' @param lwd is the line width for the estimated regression line.
#' @param col45 is the color for the 45-degree line.
#' @param col is the color for the estimated regression line.
#' @param pch is the plotting character symbol.
#' @param cex is the size of the plotting symbol.
#' @param xlim is the range of values on the x-axis.
#' @param ylim is the range of values on the y-axis.
#' @param xlab is the label on the x-axis.
#' @param ylab is the label on the y-axis.
#' @param ... additional arguments to be passed to the \code{plot}
#' function.
#' @return A plot is shown, NULL is returned
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{par}}, \code{\link{plot}}
#' @references Cunningham VJ, Rabiner EA, Slifstein M, Laruelle M (2010).
#' Measuring drug occupancy in the absence of a reference region: the Lassen
#' plot re-visited, \emph{Journal of Cerebral Blood Flow & Metababolism},
#' \bold{30}, 46-50.
#' @export plotBindingPotential
#' @importFrom graphics abline plot points
#' @importFrom stats lsfit
plotBindingPotential <- function(base, drug, lty45=2, lty=1, lwd45=2, lwd=3,
                                 col45="darkgrey", col="orange", pch=1, cex=1,
                                 xlim=range(0, base, 0.5),
                                 ylim=range(0, drug, 0.5),
                                 xlab=expression(BP[ND]^{Base}),
                                 ylab=expression(BP[ND]^{Drug}),
                                 ...) {
  plot(base, drug, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  abline(0, 1, col=col45, lty=lty45, lwd=lwd45)
  fit <- lsfit(base, drug, intercept=FALSE)
  abline(fit, lwd=lwd, col=col)
  points(base, drug, pch=pch, cex=cex)
  invisible()
}
