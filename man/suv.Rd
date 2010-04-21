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
  ...
}
\usage{
standardUptakeValue(data, mask, dose, mass)
hotSpotSUV(slice, radius=5)
totalSUV(suv, mask, z, bg, local=TRUE, mname=NULL, nt=NULL)
}
\arguments{
  \item{data}{is a multidimensional array of signal intensities.}
  \item{mask}{}
  \item{dose}{}
  \item{mass}{is the injected mass of the PET tracer (units?).}
  \item{slice}{}
  \item{radius}{is the desired hotspot radius (units?).}
  \item{suv}{is the standard uptake value (SUV).}
  \item{z}{}
  \item{bg}{}
  \item{local}{}
  \item{mname}{}
  \item{nt}{}
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

