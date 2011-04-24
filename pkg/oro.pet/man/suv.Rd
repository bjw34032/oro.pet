% Part of the oro.pet package for R
% Distributed under the BSD license: see oro.pet/COPYING
%
% $Id: $

\name{Standard Uptake Value}
\Rdversion{1.1}
\alias{standardUptakeValue}
\alias{hotSpotSUV}
\alias{totalSUV}
\title{Calculating Standard Uptake Values for PET}
\description{
  The standard uptake value (SUV) is calculated based on an 18F-FDG-PET
  acquistion using the lean body mass approximation.
}
\usage{
standardUptakeValue(data, mask, dose, mass)
hotSpotSUV(suv, radius=1, type="2D")
totalSUV(suv, mask, z, bg, local=TRUE, mname=NULL, nt=NULL)
}
\arguments{
  \item{data}{is a multidimensional array of signal intensities.}
  \item{mask}{is a multicimensional array of logical values.}
  \item{dose}{is the injected dose of the PET tracer (units?).}
  \item{mass}{is the estimated body mass of the subject (units?).}
  \item{radius}{is the desired hotspot radius (units = voxels).}
  \item{type}{is a character string (acceptable values are \code{2D} or
    \code{3D}) that determines the dimension of the hot spot (default =
    \code{2D}).}
  \item{suv}{is the standard uptake value (SUV).}
  \item{z}{is the slice index.}
  \item{bg}{is the estimated background SUV.}
  \item{local}{...}
  \item{mname}{...}
  \item{nt}{...}
}
\details{
  ...
}
\value{
  ...
}
% \references{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\seealso{
  \code{\link{leanBodyMass}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{}

