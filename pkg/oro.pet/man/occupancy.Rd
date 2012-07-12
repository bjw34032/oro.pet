\name{occupancy}
\alias{occupancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute Drug Occupancy with Approximate Standard Errors
}
\description{
  Receptor occupancy is calculated from posititron emission tomography
  (PET) data as the treatment-induced relative change in the
  concentration of available (not occupied) receptors.
}
\usage{
occupancy(base, drug, baseSE = NULL, drugSE = NULL, base.drug.corr = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base}{is the baseline binding potential (BPND).}
  \item{drug}{is the post-treatment binding potential (BPND).}
  \item{baseSE}{is the standard error for the baseline BPND.}
  \item{drugSE}{is the standard error for the post-treatment BPND.}
  \item{base.drug.corr}{is the user-specified correlation between
    baseline and post-treatment binding potentials.}
}
\details{
  Occupancy is calculated using the straightforward and well-known
  formula.  If the standard errors for the two binding potentials are
  provided, then the delta method is used to approximate the standard
  error for the estimate of occupancy.
}
\value{
  \item{OCC}{is the percent drug occupancy.}
  \item{SE}{is the approximate standard error of the parameter estimate.}
}
\references{
  Cunningham VJ, Rabiner EA, Slifstein M, Laruelle M (2010).  Measuring
  drug occupancy in the absence of a reference region: the Lassen plot
  re-visited, \emph{Journal of Cerebral Blood Flow & Metababolism},
  \bold{30}, 46-50.
  
  Passchier J, Gee A, Willemsen A, Vaalburg W, van Waarde A (2002).
  Measuring drug-related receptor occupancy with positron emission
  tomography, \emph{Methods}, \bold{27}, 278-286. 
}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[msm]{deltamethod}}
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

