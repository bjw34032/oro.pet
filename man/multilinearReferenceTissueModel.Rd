% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reference.R
\name{multilinearReferenceTissueModel}
\alias{multilinearReferenceTissueModel}
\title{The Multilinear Reference Tissue Model}
\usage{
multilinearReferenceTissueModel(
  tac,
  ref,
  time,
  tstar,
  MRTM2 = TRUE,
  k2prime = NULL
)
}
\arguments{
\item{tac}{a vector corresponding to the time activity curve from the tissue
(in Bq/mL).}

\item{ref}{a vector corresponding to the time activity curve from the
reference region (in Bq/mL).}

\item{time}{a vector of average frame times (in minutes).}

\item{tstar}{the time (in minutes) where the linear relationship between the
response and covariates may be assumed to be true.}

\item{MRTM2}{a logical value that selects the three-parameter model (MRTM)
or the two-parameter model (MRTM2), where k2prime is fixed.}

\item{k2prime}{the value of k2prime that has been fixed.}
}
\value{
\item{BP}{Binding potential} \item{BP.error}{Approximate standard
error of the binding potential} \item{R1}{Ratio of the volumes of
distrubution for the tissue and reference region (assumes a one-tissue model
is valid)} \item{R1.error}{Approximate standard error for the ratio}
\item{k2}{Clearance rate constant from the tissue to plasma (assumes a
one-tissue model is valid)} \item{k2.error}{Approximate standard error for
k2} \item{X}{Design matrix used in the linear regression}
\item{beta}{Regression coefficients}
}
\description{
The multilinear reference tissue model (MRTM) estimates the binding
potential from an observed time activity curve without the need for arterial
sampling.  Instead, a second time activity curve must be provided from a
suitable reference region where there is negligible binding.
}
\details{
See the references.

The numeric integration required to construct the design matrix is performed
by interpolating the time activity curves, both for the tissue and reference
region, to one-second resolution and then performing the \code{cumsum}
operation on them.

Given the nonlinear relationship between binding potential and the
regression parameters, the \code{deltamethod} is used to approximate its
standard error.
}
\references{
Ichise, M., Ballinger, J.R., Golan, H., Vines, D., Luong, A.,
Tsai, S. and Kung, H.F. (1996) Noninvasive quantification of dopamine D2
receptors with iodine-123-IBF SPECT, \emph{Journal of Nuclear Medicine},
\bold{37}, 513-520.

Ichise, M., Liow, J.-S., Lu, J.-Q., Takano, A., Model, K., Toyama, H.,
Suhara, T., Suzuki, K., Innis, R.B., Carson, R.E. (2003) Linearized
reference tissue parametric imaging methods: Application to [11C]DASB
positron emission tomography studies of the serotonin transporter in human
brain, \emph{Journal of Cerebral Blood Flow & Metabolism}, \bold{23},
1096-1112.
}
\seealso{
\code{\link{cumsum}}, \code{\link[msm]{deltamethod}}
}
\author{
Brandon Whitcher \email{bwhitcher@gmail.com}
}
