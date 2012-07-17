\name{plotBindingPotential}
\alias{plotBindingPotential}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Baseline Versus Post-Treatment Binding Potentials
}
\description{
  Inspired by the Lassen plot (Cunningham et al., 2010) this is a
  straightforward graphical summary of pre-treatment versus
  post-treatment binding potentials for a single subject across multiple
  brain regions.
}
\usage{
plotBindingPotential(base, drug, lty45=2, lty=1, lwd45=2, lwd=3,
                     col45="darkgrey", col="orange", pch=1, cex=1,
                     xlim=range(0, base, 0.5), ylim=range(0, drug, 0.5),
                     xlab=expression(BP[ND]^{Base}),
                     ylab=expression(BP[ND]^{Drug}), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base}{is the vector of baseline binding potentials across brain
    regions.}
  \item{drug}{is the vector of post-treatment binding potentials across
    brain regions.}
  \item{lty45}{is the line type for the 45-degree line.}
  \item{lty}{is the line type for the estimated regression line.}
  \item{lwd45}{is the line width for the 45-degree line.}
  \item{lwd}{is the line width for the estimated regression line.}
  \item{col45}{is the color for the 45-degree line.}
  \item{col}{is the color for the estimated regression line.}
  \item{pch}{is the plotting character symbol.}
  \item{cex}{is the size of the plotting symbol.}
  \item{xlim}{is the range of values on the x-axis.}
  \item{ylim}{is the range of values on the y-axis.}
  \item{xlab}{is the label on the x-axis.}
  \item{ylab}{is the label on the y-axis.}
  \item{\ldots}{for additional arguments to be passed to the \code{plot}
    fucntion.}
}
\details{
  See the reference below.
}
\value{
  A plot.
}
\references{
  Cunningham VJ, Rabiner EA, Slifstein M, Laruelle M (2010).  Measuring
  drug occupancy in the absence of a reference region: the Lassen plot
  re-visited, \emph{Journal of Cerebral Blood Flow & Metababolism},
  \bold{30}, 46-50.
}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{par}}, \code{\link{plot}}
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

