\name{expConv}
\alias{expConv}
\title{
  Empirical Convolution Between an Input Function and a Single Exponential
}
\description{
  Computationally efficient method to convolve a vector of observations
  and a single exponential function with two parameters.
}
\usage{
expConv(input, k1, k2)
}
\arguments{
  \item{input}{is the so-called input function.}
  \item{k1}{is the scaling parameter in the single exponential function.}
  \item{k2}{is the decay parameters in the single exponential function.}
}
\details{
  Assuming the input function has been sampled (or interpolated) to a
  high temporal resolutions, say one Hertz, a simple for loop is used to
  perform the convolution.  
}
\value{
  The vector containing the result from the convolution operation.
}
% \references{}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
% \seealso{}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
