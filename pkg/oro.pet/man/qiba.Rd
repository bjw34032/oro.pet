% Part of the oro.pet package for R
% Distributed under the BSD license: see oro.pet/COPYING
%
% $Id: $

\name{QIBA}
\Rdversion{1.1}
\alias{standardUptakeValue}
\alias{standardUptakeValue,array-method}
\alias{activityConcentration}
\alias{activityConcentration,array-method}
\title{Calculating SUVs for PET Using QIBA Pseudocode}
\description{
  The standard uptake value (SUV) is calculated based on an 18F-FDG-PET
  acquistion using ancillary information contained in the DICOM data. 
}
\usage{
\S4method{standardUptakeValue}{array}(pixelData, mask = NULL, CSV = NULL,
                    seriesNumber = NULL, method = c("qiba", "user"),
                    prior = NULL, decayedDose = NULL)
\S4method{activityConcentration}{array}(pixelData, CSV = NULL, seriesNumber = NULL,
                      method = "qiba")
}
\arguments{
  \item{pixelData}{is a multidimensional array of signal intensities of
    class \code{nifti}.}
  \item{mask}{is a multidimensional array of logical values (only used
    when \code{method = "user"}).}
  \item{CSV}{is a \code{data.frame} that is the output from
    \code{dicom.table} and contains all necessary DICOM header fields.}
  \item{seriesNumber}{is the SeriesNumber that corresponds to the PET
    acquisition.}
  \item{method}{takes on two possible values (\code{qiba} and
    \code{user}), where QIBA pseudocode is used to calculate the SUVs or
    user-defined parameters are used.}
  \item{prior}{is a list of DICOM header field names that are necessary
    for the SUV calculation under \code{method = "user"} or may be used
    to replace values from the DICOM header information when
    \code{method = "qiba"}.}
  \item{decayedDose}{is the amount of the RadionuclideTotalDose after
    being corrected for residual dose in the syringe.  This value is NOT
    usually corrected in the DICOM data.}
}
\details{
  \bold{[MORE]}

  Note, for GE scanners it is common for the RescaleSlope DICOM field to
  vary on a slice-by-slice basis.  This is taken into account if a GE
  scanner is detected from the Modality DICOM field.  However, the
  InstanceNumber is used to reorder the slices so they match the
  incoming NIfTI file of PixelData.  If this is not correct it may be
  necessary to manually re-order the RescaleSlope field in the CSV data
  frame so that the activity concentration is calculated correctly.
}
\value{
  A list containing the following items
  \item{SUVbw}{is a multidimensional array, the same dimension as
    \code{pixelData}, that contains the standard uptake values.}
  \item{hdr}{is a list of DICOM header fields used in the SUV
  calculation.}
  \item{decayTime}{is the decay time calculated from the DICOM header
    information.}
  \item{decayedDose}{is the RadionuclideTotalDose, if taken from the
    DICOM header information, or the user-specified value.}
  \item{SUVbwScaleFactor}{is
    \eqn{\mbox{PatientsWeight}\cdot1000/\mbox{decayedDose}}.}
}
\references{
  \url{http://qibawiki.rsna.org/index.php?title=Standardized_Uptake_Value_(SUV)}
}
\author{Brandon Whitcher \email{bwhitcher@gmail.com}}
\seealso{
  \code{\link[oro.dicom]{dicomTable}}, \code{\link[oro.nifti]{nifti}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{}

