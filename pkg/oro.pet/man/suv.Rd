% Part of the oro.pet package for R
% Distributed under the BSD license: see oro.pet/COPYING
%
% $Id: $

\name{Summarizing SUVs}
\Rdversion{1.1}
%\alias{standardUptakeValue}
\alias{hotSpotSUV}
\alias{totalSUV}
\title{Summarizing SUVs for PET}
\description{
  The standard uptake value (SUV) is summarized using the hotspot method
  or by calculating total volume of the high values.
}
\usage{
%standardUptakeValue(data, mask, dose, mass)
hotSpotSUV(suv, radius = 10, type = "3D")
totalSUV(suv, mask, z, bg, local = TRUE)
}
\arguments{
  \item{suv}{is the standard uptake value (SUV).}
  \item{radius}{is the desired hotspot radius (units = voxels).}
  \item{type}{is a character string (acceptable values are \code{2D} or
    \code{3D}) that determines the dimension of the hot spot (default =
    \code{3D}).}
  \item{mask}{is a multidimensional array of logical values.}
  \item{z}{is the slice index.}
  \item{bg}{is the estimated background SUV.}
  \item{local}{is a logical value.}
}
\details{
  ...
}
\value{
  ...
}
% \references{}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
\seealso{
  \code{\link{leanBodyMass}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{}

