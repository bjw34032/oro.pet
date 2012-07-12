% Part of the oro.pet package for R
% Distributed under the BSD license: see oro.pet/COPYING
%
% $Id: $

\name{LeanBodyMass}
\Rdversion{1.1}
\alias{leanBodyMass}
\title{Calculating the Lean Body Mass}
\description{
  The lean body mass (LBM) is calculated according to the formula
  \deqn{1.1\cdot\mbox{weight}-128\cdot(\mbox{weight}/\mbox{height})^2}
  if male and
  \deqn{1.07\cdot\mbox{weight}-148\cdot(\mbox{weight}/\mbox{height})^2}
  if female.
}
\usage{
leanBodyMass(height, weight, gender)
}
\arguments{
  \item{height}{is a vector of heights in centimeters.}
  \item{weight}{is a vector of weights in kilograms.}
  \item{gender}{is a character vector (may be of length one) with the
    value \dQuote{male} or \dQuote{female}.}
}
%\details{}
\value{
  Vector of lean body mass values in kilograms.
}
\references{
  Sugawara, Y., K. R. Zasadny, A. W. Neuhoff, R. L. Wahl (1999)
  Reevaluation of the Standardized Uptake Value for FDG: Variations with
  Body Weight and Methods for Correction, \emph{Radiology} \bold{213}:
  521--525.
}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
\seealso{
  \code{\link{standardUptakeValue}}
}
\examples{
n <- 11
h <- seq(200, 150, length=n)
w <- seq(80, 120, length=n)
cbind(h, w, leanBodyMass(h, w, "male"), leanBodyMass(h, w, "female"))
}
% \keyword{ ~kwd1 }

