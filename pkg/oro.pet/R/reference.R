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

simplifiedReferenceTissueModel <- function(tac, ref, time, SRTM2=TRUE,
                                           k2prime=NULL,
                                           guess=c("R1"=0.5, "k2"=0.01),
                                           control=nls.lm.control()) {
  require("minpack.lm")
  require("msm")
  func <- function(theta, signal, time, ref, k2prime) {
    out <- signal - extended.empirical(time, theta, ref, k2prime)
    out[!is.na(out)]
  }
  out <- nls.lm(par=guess, fn=func, control=control, signal=tac, time=time,
                ref=ref, k2prime=k2prime)
  srtm2 <- list(BP = out$par[1] * k2prime / out$par[2] - 1,
                R1 = out$par[1],
                k2 = out$par[2])
  srtm2 <- lapply(srtm2, as.numeric)
  rdf <- length(out$fvec) - length(coef(out))
  varcovmat <- (out$deviance / rdf) * chol2inv(chol(out$hessian))
  c(srtm2, BP.error = deltamethod(~ x1 * k2prime / x2, out$par, varcovmat),
    R1.error = sqrt(varcovmat[1,1]), k2.error = sqrt(varcovmat[2,2]))
}

multilinearReferenceTissueModel <- function(tac, ref, time, tstar,
                                            MRTM2=TRUE, k2prime=NULL) {
  require("msm")
  ## Numeric integration
  time.in.sec <- seq(min(0, time * 60), ceiling(max(time * 60)), by=1)
  sec <- list(tac = approx(c(0, time) * 60, c(0, tac), time.in.sec)$y,
              ref = approx(c(0, time) * 60, c(0, ref), time.in.sec)$y)
  ## Fit linear model and estimate parameters
  X <- cbind(cumsum(sec$ref) + sec$ref / (k2prime / 60), # in seconds
             cumsum(sec$tac))[time * 60,]
  dimnames(X) <- list(NULL, paste("gamma", 1:2, sep=""))
  index <- time > tstar
  fit <- lsfit(X[index,], tac[index], intercept=FALSE)
  gamma <- fit$coefficients
  mrtm2 <- list(BP = - (gamma[1] / gamma[2] + 1),
                R1 = gamma[1] / k2prime,
                k2 = - gamma[2])
  mrtm2 <- lapply(mrtm2, as.numeric)
  ## Construct variance-covariance matrix for regression parameters
  rdf <- length(residuals(fit)) - length(coefficients(fit))
  varcovbeta <- (sum(residuals(fit)^2) / rdf) * chol2inv(chol(t(X) %*% X))
  ## Construct list output with approximate standard errors
  list(BP = as.numeric(mrtm2$BP),
       BP.error = deltamethod(~ x1 / x2, gamma, varcovbeta),
       R1 = as.numeric(mrtm2$R1), R1.error = as.numeric(sqrt(varcovbeta[1,1])),
       k2 = as.numeric(mrtm2$k2), k2.error = as.numeric(sqrt(varcovbeta[2,2])),
       X = X, beta = gamma)
}
