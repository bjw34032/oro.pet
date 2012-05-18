\name{compartmentalModel}
\alias{compartmentalModel}
\title{
  Compartmental Models for Kinetic Parameter Estimation
}
\description{
  A selection of parametric models are provided that combine a
  compartmental model for tissue and empirical versions of the arterial
  input function or reference region time activity curve.
}
\usage{
compartmentalModel(type)
}
\arguments{
  \item{type}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \describe{
      \item{\dQuote{srtm}}{Simplified Reference Tissue Model}
    }
  }
}
\details{
  Parametric models from the PET literature are provided to the user for
  kinetic parameter estimation.
}
\value{
  A function.
}
\references{

}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
\seealso{
  \code{\link{simplifiedReferenceTissueModel}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
