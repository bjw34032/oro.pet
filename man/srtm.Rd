\name{simplifiedReferenceTissueModel}
\alias{simplifiedReferenceTissueModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The Simplified Reference Tissue Model
}
\description{
  The simplified reference tissue model (SRTM) estimates the binding
  potential from an observed time activity curve without the need for
  aterial sampling.  It assumes a one-tissue compartment model to
  describe the influx and efflux in the tissue region of interest and
  the reference region.
}
\usage{
simplifiedReferenceTissueModel(tac, ref, time, SRTM2 = TRUE, k2prime = NULL,
                               guess = c(R1 = 0.5, k2 = 0.01),
                               control = minpack.lm::nls.lm.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tac}{a vector corresponding to the time activity curve from the
    tissue (in Bq/mL).}
  \item{ref}{a vector corresponding to the time activity curve from the
    reference region (in Bq/mL).}
  \item{time}{a vector of average frame times (in minutes).}
  \item{SRTM2}{a logical value that selects the three-parameter model
    (SRTM) or the two-parameter model (SRTM2), where k2prime is fixed.}
  \item{k2prime}{the value of k2prime that has been fixed.}
  \item{guess}{values for the inital parameter estimates for R1 and k2.}
  \item{control}{a list of parameters used by \code{nls.lm.control}
    that are set by default, but may be customized by the user.}
}
\details{
  See the references.

  The model has been parameterized in the manner of Wu and Carson
  (2002).  That is, the nonlinear regression estimates R1, k2 and k'2
  for the three-parameter model (SRTM) and R1 and k2 for the
  two-parameter model (SRTM2).

  The convolution is performed after interpolating the time activity
  curves, both for the tissue and the reference region, to one-second
  resolution then downsampling them back to the original sampling rate.
}
\value{
  \item{BP}{Binding potential}
  \item{R1}{Ratio of the volumes of distrubution for the tissue and
    reference region}
  \item{k2}{Clearance rate constant from the tissue to plasma}
  \item{BP.error}{Approximate standard error of the binding potential}
  \item{R1.error}{Approximate standard error for the ratio}
  \item{k2.error}{Approximate standard error for k2}
}
\references{
  Lammertsma, A.A. and Hume, S.P. (1996) Simplified reference tissue
  model for PET receptor studies, \emph{NeuroImage}, \bold{4}, 153-158.

  Wu, Y. and Carson, R.E. (2002) Noise reduction in the simplified
  reference tissue model for neuroreceptor functional imaging,
  \emph{Journal of Cerebral Blood Flow \& Metabolism}, \bold{22},
  1440-1452.
}
\author{Brandon Whitcher \email{b.whitcher@gmail.com}}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[msm]{deltamethod}}, \code{\link{expConv}},
  \code{\link[minpack.lm]{nls.lm}}
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

