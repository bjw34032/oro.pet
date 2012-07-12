\name{hillEquation}
\alias{hillEquation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimation of the Half Maximal Inhibitory Concentration
}
\description{
  The half maximal inhibitory concentration (IC50) is a measure of the
  effectiveness of a compound in inhibiting biological or biochemical
  function.  This quantitative measure indicates how much of a
  particular drug or other substance (inhibitor) is needed to inhibit a
  given biological process (or component of a process) by half.
}
\usage{
hillEquation(conc, occ, guess = c(1, 100), control = nls.lm.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conc}{a vector of drug concentrations in plasma (example units
    are ng/mL).}
  \item{occ}{a vector of PET occupancy values that correspond to the
    measured drug concentrations in plasma.}
  \item{guess}{a length-two vector of starting values for the nonlinear
    optimization.}
  \item{control}{is a list of parameters used by \code{nls.lm.control}
    that are set by default, but may be customized by the user.}
}
\details{
  See reference(s).

  In this version of the function the maximal occupancy (rmax) is
  estimated automatically.  This should be optional.
}
\value{
  \item{IC50}{Half maximal inhibitory concentration}
  \item{rmax}{Estimated maximal occupancy}
  \item{IC50SE}{Approximate standard error for IC50}
  \item{rmaxSE}{Approximate standard erorr for rmax}
  \item{hessian}{Hessian matrix from the Levenburg-Marquardt procedure}
  \item{info}{Return value from the Levenburg-Marquardt procedure}
  \item{deviance}{Deviance from the Levenburg-Marquardt procedure}
  \item{message}{Text message from the Levenburg-Marquardt procedure}
}
\references{
  \href{Hill Equation}{http://en.wikipedia.org/wiki/Hill_equation_(biochemistry)}
  \href{IC50}{http://en.wikipedia.org/wiki/IC50}
}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[minpack.lm]{nls.lm}}
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

