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



#' The Simplified Reference Tissue Model
#' 
#' The simplified reference tissue model (SRTM) estimates the binding potential
#' from an observed time activity curve without the need for aterial sampling.
#' It assumes a one-tissue compartment model to describe the influx and efflux
#' in the tissue region of interest and the reference region.
#' 
#' See the references.
#' 
#' The model has been parameterized in the manner of Wu and Carson (2002).
#' That is, the nonlinear regression estimates R1, k2 and k'2 for the
#' three-parameter model (SRTM) and R1 and k2 for the two-parameter model
#' (SRTM2).
#' 
#' The convolution is performed after interpolating the time activity curves,
#' both for the tissue and the reference region, to one-second resolution then
#' downsampling them back to the original sampling rate.
#' 
#' @param tac a vector corresponding to the time activity curve from the tissue
#' (in Bq/mL).
#' @param ref a vector corresponding to the time activity curve from the
#' reference region (in Bq/mL).
#' @param time a vector of average frame times (in minutes).
#' @param SRTM2 a logical value that selects the three-parameter model (SRTM)
#' or the two-parameter model (SRTM2), where k2prime is fixed.
#' @param k2prime the value of k2prime that has been fixed.
#' @param guess values for the inital parameter estimates for R1 and k2.
#' @param control a list of parameters used by \code{nls.lm.control} that are
#' set by default, but may be customized by the user.
#' @return \item{BP}{Binding potential} \item{R1}{Ratio of the volumes of
#' distrubution for the tissue and reference region} \item{k2}{Clearance rate
#' constant from the tissue to plasma} \item{BP.error}{Approximate standard
#' error of the binding potential} \item{R1.error}{Approximate standard error
#' for the ratio} \item{k2.error}{Approximate standard error for k2}
#' @author Brandon Whitcher \email{b.whitcher@@gmail.com}
#' @seealso \code{\link[msm]{deltamethod}}, \code{\link{expConv}},
#' \code{\link[minpack.lm]{nls.lm}}
#' @references Lammertsma, A.A. and Hume, S.P. (1996) Simplified reference
#' tissue model for PET receptor studies, \emph{NeuroImage}, \bold{4}, 153-158.
#' 
#' Wu, Y. and Carson, R.E. (2002) Noise reduction in the simplified reference
#' tissue model for neuroreceptor functional imaging, \emph{Journal of Cerebral
#' Blood Flow & Metabolism}, \bold{22}, 1440-1452.
#' @importFrom minpack.lm nls.lm nls.lm.control
simplifiedReferenceTissueModel <- function(tac, ref, time, SRTM2=TRUE,
                                           k2prime=NULL,
                                           guess=c("R1"=0.5, "k2"=0.01),
                                           control=minpack.lm::nls.lm.control()) {
    # require("msm")
  func.model <- compartmentalModel(ifelse(SRTM2, "srtm2", "srtm"))
  func <- function(theta, signal, time, ref, k2prime) {
    vec <- signal - func.model(time, theta, ref, k2prime)
    vec[!is.na(vec)]
  }
  nlls <- minpack.lm::nls.lm(par=guess, fn=func, control=control, signal=tac, 
                             time=time, ref=ref, k2prime=k2prime)
  ## Construct variance-covariance matrix for regression parameters
  rdf <- length(nlls$fvec) - length(coef(nlls))
  varcovmat <- (nlls$deviance / rdf) * chol2inv(chol(nlls$hessian))
  ## Construct list output with approximate standard errors
  list(BP = as.numeric(nlls$par[1] * k2prime / nlls$par[2] - 1),
       R1 = as.numeric(nlls$par[1]),
       k2 = as.numeric(nlls$par[2]),
       BP.error = msm::deltamethod(~ x1 * k2prime / x2, nlls$par, varcovmat),
       R1.error = sqrt(varcovmat[1,1]),
       k2.error = sqrt(varcovmat[2,2]),
       hessian = nlls$hessian, info = nlls$info, deviance = nlls$deviance,
       message = nlls$message)
}



#' The Multilinear Reference Tissue Model
#' 
#' The multilinear reference tissue model (MRTM) estimates the binding
#' potential from an observed time activity curve without the need for arterial
#' sampling.  Instead, a second time activity curve must be provided from a
#' suitable reference region where there is negligible binding.
#' 
#' See the references.
#' 
#' The numeric integration required to construct the design matrix is performed
#' by interpolating the time activity curves, both for the tissue and reference
#' region, to one-second resolution and then performing the \code{cumsum}
#' operation on them.
#' 
#' Given the nonlinear relationship between binding potential and the
#' regression parameters, the \code{deltamethod} is used to approximate its
#' standard error.
#' 
#' @param tac a vector corresponding to the time activity curve from the tissue
#' (in Bq/mL).
#' @param ref a vector corresponding to the time activity curve from the
#' reference region (in Bq/mL).
#' @param time a vector of average frame times (in minutes).
#' @param tstar the time (in minutes) where the linear relationship between the
#' response and covariates may be assumed to be true.
#' @param MRTM2 a logical value that selects the three-parameter model (MRTM)
#' or the two-parameter model (MRTM2), where k2prime is fixed.
#' @param k2prime the value of k2prime that has been fixed.
#' @return \item{BP}{Binding potential} \item{BP.error}{Approximate standard
#' error of the binding potential} \item{R1}{Ratio of the volumes of
#' distrubution for the tissue and reference region (assumes a one-tissue model
#' is valid)} \item{R1.error}{Approximate standard error for the ratio}
#' \item{k2}{Clearance rate constant from the tissue to plasma (assumes a
#' one-tissue model is valid)} \item{k2.error}{Approximate standard error for
#' k2} \item{X}{Design matrix used in the linear regression}
#' \item{beta}{Regression coefficients}
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{cumsum}}, \code{\link[msm]{deltamethod}}
#' @references Ichise, M., Ballinger, J.R., Golan, H., Vines, D., Luong, A.,
#' Tsai, S. and Kung, H.F. (1996) Noninvasive quantification of dopamine D2
#' receptors with iodine-123-IBF SPECT, \emph{Journal of Nuclear Medicine},
#' \bold{37}, 513-520.
#' 
#' Ichise, M., Liow, J.-S., Lu, J.-Q., Takano, A., Model, K., Toyama, H.,
#' Suhara, T., Suzuki, K., Innis, R.B., Carson, R.E. (2003) Linearized
#' reference tissue parametric imaging methods: Application to [11C]DASB
#' positron emission tomography studies of the serotonin transporter in human
#' brain, \emph{Journal of Cerebral Blood Flow & Metabolism}, \bold{23},
#' 1096-1112.
#' @importFrom stats coefficients residuals
#' @importFrom msm deltamethod
multilinearReferenceTissueModel <- function(tac, ref, time, tstar,
                                            MRTM2=TRUE, k2prime=NULL) {
  # require("msm")
  ## Numeric integration
  time.in.sec <- seq(min(0, time * 60), ceiling(max(time * 60)), by=1)
  sec <- list(tac = approx(c(0, time) * 60, c(0, tac), time.in.sec)$y,
              ref = approx(c(0, time) * 60, c(0, ref), time.in.sec)$y)
  X <- cbind(cumsum(sec$ref) + sec$ref / (k2prime / 60), # in seconds
             cumsum(sec$tac))[time * 60,]
  dimnames(X) <- list(NULL, paste("gamma", 1:2, sep=""))
  ## Fit linear model and estimate parameters
  index <- time > tstar
  fit <- lsfit(X[index,], tac[index], intercept=FALSE)
  gamma <- fit$coefficients
  ## Construct variance-covariance matrix for regression parameters
  rdf <- length(residuals(fit)) - length(coefficients(fit))
  varcovbeta <- (sum(residuals(fit)^2) / rdf) * chol2inv(chol(t(X) %*% X))
  ## Construct list output with approximate standard errors
  list(BP = as.numeric(- (gamma[1] / gamma[2] + 1)),
       R1 = as.numeric(gamma[1] / k2prime),
       k2 = as.numeric(- gamma[2]),
       BP.error = msm::deltamethod(~ x1 / x2, gamma, varcovbeta),
       R1.error = as.numeric(sqrt(varcovbeta[1,1])),
       k2.error = as.numeric(sqrt(varcovbeta[2,2])),
       X = X, beta = gamma)
}
