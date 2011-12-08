##
## Copyright (c) 2011, Brandon Whitcher
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

#############################################################################
## setGeneric("activityConcentration")
#############################################################################

setGeneric("activityConcentration",
           function(pixelData,  ...) standardGeneric("activityConcentration"))
setMethod("activityConcentration", signature(pixelData = "array"), 
	  function(pixelData, CSV = NULL, seriesNumber = NULL, method = "qiba")
          .petWrapper("activityConcentration", pixelData, CSV,
                      seriesNumber, method))

.activityConcentration <- function(pixelData, CSV=NULL, seriesNumber=NULL,
                                   method="qiba") {
  if (is.null(CSV)) {
    stop("CSV file from dicom.table() is required")
  }
  if (is.null(seriesNumber)) {
    stop("SeriesNumber of PET acquisition is required")
  }
  ##
  ## rescale FDG-PET using Y = mX + b
  ##
  csv <- CSV[CSV[, "X0020.0011.SeriesNumber"] == seriesNumber, ]
  acquisitionTime <- csv[, "X0008.0032.AcquisitionTime"]
  instanceNumber <- csv[, "X0020.0013.InstanceNumber"]
  instanceUID <- csv[, "X0008.0018.SOPInstanceUID"]
  if (length(instanceNumber) != length(unique(instanceNumber))) {
    warning("InstanceNumber is not a unique identifier")
  }
  ino <- order(instanceNumber) # the order of instanceNumber
  rescaleIntercept <- csv[, "X0028.1052.RescaleIntercept"] # b
  rescaleSlope <- csv[, "X0028.1053.RescaleSlope"] # m
  rI <- rS <- matrix(0, nsli(pixelData), ntim(pixelData))
  activity <- array(0, dim(pixelData))
  nslices <- nsli(pixelData)
  for (i in 1:length(rescaleSlope)) {
    z <- (i - 1) %% nslices + 1
    w <- (i - 1) %/% nslices + 1
    j <- ino[i] # each GE correction is slice-specific!
    rI[z,w] <- rescaleIntercept[j]
    rS[z,w] <- rescaleSlope[j]
    activity[,,z,w] <- pixelData[,,z,w] * rescaleSlope[j] + rescaleIntercept[j]
  }
  list(activity = activity, rescaleIntercept = rescaleIntercept,
       rescaleSlope = rescaleSlope, resclaeInterceptMatrix = rI,
       rescaleSlopeMatrix = rS)
}

#############################################################################
## setGeneric("standardUptakeValue")
#############################################################################

setGeneric("standardUptakeValue",
           function(pixelData,  ...) standardGeneric("standardUptakeValue"))
setMethod("standardUptakeValue", signature(pixelData = "array"), 
	  function(pixelData, mask = NULL, CSV = NULL,
                   seriesNumber = NULL, method = c("qiba", "user"),
                   prior = NULL, decayedDose = NULL)
          .petWrapper("standardUptakeValue", pixelData, mask, CSV,
                      seriesNumber, method, prior, decayedDose))

.standardUptakeValue <- function(pixelData, mask=NULL, CSV=NULL,
                                 seriesNumber=NULL, method=c("qiba","user"),
                                 prior=NULL, decayedDose=NULL) {
  ## user = quick-and-dirty method
  ## qiba = QIBA-based pseudocode that requires DICOM data
  ## 
  if (method[1] == "user") {
    suv <- array(NA, dim(pixelData))
    suv[mask] <- pixelData[mask] / prior$dose * prior$mass
    return(suv)
  } else {
    if (is.null(CSV)) {
      stop("CSV file from dicom.table() is required")
    }
    if (is.null(seriesNumber)) {
      stop("SeriesNumber of PET acquisition is required")
    }
    ## SUV cannot be calculated if any of the specified DICOM attributes
    ## are missing or empty or zero
    
    ## Obtain DICOM header information from .csv file
    colSerNum <- "X0020.0011.SeriesNumber"
    hdr <- vector(mode="list")
    hdr$correctedImage <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                     "X0028.0051.CorrectedImage"])
    hdr$decayCorrection <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                      "X0054.1102.DecayCorrection"])
    hdr$units <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                            "X0054.1001.Units"])
    hdr$seriesDate <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                 "X0008.0021.SeriesDate"])
    hdr$seriesDate <- as.character(hdr$seriesDate)
    hdr$seriesTime <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                 "X0008.0031.SeriesTime"])
    hdr$seriesTime <- as.character(hdr$seriesTime)
    hdr$acquisitionDate <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                      "X0008.0022.AcquisitionDate"])
    hdr$acquisitionDate <- as.character(hdr$acquisitionDate)
    hdr$acquisitionTime <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                      "X0008.0032.AcquisitionTime"])
    hdr$acquisitionTime <- as.character(hdr$acquisitionTime)
    ## RadiopharmaceuticalStartTime is in RadiopharmaceuticalInformationSequence
    hdr$radiopharmaceuticalStartTime <-
      unique(CSV[CSV[,colSerNum] == seriesNumber,
                 "X0054.0016.0018.1072.RadiopharmaceuticalStartTime"])
    ## RadionuclideTotalDose is in RadiopharmaceuticalInformationSequence
    hdr$radionuclideTotalDose <-
      unique(CSV[CSV[,colSerNum] == seriesNumber,
                 "X0054.0016.0018.1074.RadionuclideTotalDose"])
    ## RadionuclideHalfLife is in RadiopharmaceuticalInformationSequence
    hdr$radionuclideHalfLife <-
      unique(CSV[CSV[,colSerNum] == seriesNumber,
                 "X0054.0016.0018.1075.RadionuclideHalfLife"])
    hdr$gePrivateScanDateAndTime <-
      unique(CSV[CSV[,colSerNum] == seriesNumber,
                 "X0009.100D.Missing"]) # "GEMS_PETD_01"
    hdr$patientsWeight <- unique(CSV[CSV[,colSerNum] == seriesNumber,
                                     "X0010.1030.PatientsWeight"])
    hdr$rescaleIntercept <- CSV[CSV[,colSerNum] == seriesNumber,
                                "X0028.1052.RescaleIntercept"]
    hdr$rescaleSlope <- CSV[CSV[,colSerNum] == seriesNumber,
                            "X0028.1053.RescaleSlope"]
    hdr$instanceNumber <- CSV[CSV[,colSerNum] == seriesNumber,
                              "X0020.0013.InstanceNumber"]
    if (length(hdr$instanceNumber) != length(unique(hdr$instanceNumber))) {
      warning("InstanceNumber is not a unique identifier")
    }
    ino <- order(hdr$instanceNumber) # the order of instanceNumber
    ## Replace "hdr" information with user input
    if (! is.null(prior)) {
      for (i in 1:length(prior)) {
        ## must double-check stuff!
        j <- which(names(hdr) %in% names(prior)[i])
        hdr[[j]] <- prior[[i]]
      }
    }
    ## BEGIN pseudocode
    if (grepl("ATTN", hdr$correctedImage) &&
        grepl("DEC*Y", hdr$correctedImage) &&
        hdr$decayCorrection == "START") {
      if (hdr$units =="BQML") {
        if (str2date(hdr$seriesDate) <= str2date(hdr$acquisitionDate) &&
            all(str2time(hdr$seriesTime)$time <= str2time(hdr$acquisitionTime)$time))  {
          scanDate <- hdr$seriesDate
          scanTime <- hdr$seriesTime
        } else { # may be post-processed series
          stop("GE private scan Date and Time")
          ## scan Date and Time = GE private scan Date and Time
          ## (0x0009,0x100d,"GEMS_PETD_01")
        }
        startTime <- hdr$radiopharmaceuticalStartTime
        ## start Date is not explicit ... assume same as Series Date; but
        ## consider spanning midnight
        decayTime <- str2time(scanTime)$time - str2time(startTime)$time # seconds
        halfLife <- hdr$radionuclideHalfLife  # seconds
        ## Radionuclide Total Dose is NOT corrected for residual dose in
        ## syringe, which is ignored here ...
        injectedDose <- hdr$radionuclideTotalDose # Bq
        if (is.null(decayedDose)) {
          decayedDose <- injectedDose * 2^(-decayTime / halfLife)
        }
        weight <- hdr$patientsWeight # in kg
        SUVbwScaleFactor <- (weight * 1000 / decayedDose)
      } else {
        if (hdr$units == "CNTS") {
          stop("Philips private scale factor")
          ## SUVbwScaleFactor = Philips private scale factor
          ## (0x7053,0x1000,"Philips PET Private Group")
          ##
          ## if (0x7053,0x1000) not present, but (0x7053,0x1009) is present,
          ## then (0x7053,0x1009) * Rescale Slope
          ## scales pixels to Bq/ml, and proceed as if Units are BQML
        } else {
          if (hdr$units == "GML") {
            ## assumes that GML indicates SUVbw instead of SUVlbm
            SUVbwScaleFactor <- 1.0  
          }
        }
      }
      ## Rescale Intercept is required to be 0 for PET, but use it just in case
      ##
      ## Rescale slope may vary per slice (GE), and cannot be assumed to
      ## be constant for the entire volume
      ##
      ## I have reordered RescaleIntercept and RescaleSlope by InstanceNumber
      ## (may only be valid for GE scanners)
      SUVbw <- array(0, dim(pixelData))
      nslices <- nsli(pixelData)
      for (i in 1:length(hdr$rescaleSlope)) {
        z <- (i - 1) %% nslices + 1
        w <- (i - 1) %/% nslices + 1
        j <- ino[i] # each GE correction is slice-specific!
        SUVbw[,,z,w] <- (pixelData[,,z,w] * hdr$rescaleSlope[j] * SUVbwScaleFactor) + hdr$rescaleIntercept[j]
      }
      list(SUVbw=SUVbw, hdr=hdr, decayTime=decayTime,
           decayedDose=decayedDose, SUVbwScaleFactor=SUVbwScaleFactor)
    } else {
      stop("ATTN, START, DEC*Y... failed")
    }
    ## END pseudocode
  }
}


