LeanBodyMass <- function(height, weight, gender="male") {
  ## weight in kg
  ## height in cm
  switch(gender,
         "male" =  1.1 * weight - 128 * (weight / height)^2,
         "female" = 1.07 * weight - 148 * (weight / height)^2,
         stop("Incorrect gender specified"))
}

StandardUptakeValue <- function(data, mask, dose, mass) {
  ## Make sure the units cancel!
  ifelse(mask, data / dose * mass, NA)
}

HotSpotSUV <- function(slice, radius) {
  circle <- function(X, Y, c1, c2, r=1) {
    x <- matrix(1:X, X, Y, byrow=TRUE)
    y <- matrix(1:Y, X, Y)
    (c1 - x)^2 + (c2 - y)^2 <= r^2
  }
  m <- which(slice == max(slice, na.rm=TRUE), arr.ind=TRUE)
  slice.max <- slice[m]
  hs.mask <- circle(ncol(slice), nrow(slice), m[2], m[1], radius)
  hs.mask <- ifelse(hs.mask > 0, TRUE, FALSE)
  mean(slice[hs.mask], na.rm=TRUE)
}

TotalSUV <- function(suv, mask, z, bg, local=TRUE, mname=NULL, nt=NULL) {
  j <- which(suv == max(suv, na.rm=TRUE), arr.ind=TRUE)
  suv.max <- suv[j]
  Z <- length(z)
  vol <- numeric(Z)
  vol.mask <- array(FALSE, dim(suv))
  suv.total <- vector("list", Z)
  for (i in 1:Z) {
    sli <- suv[, , z[i]]
    if (local) {
      k <- which(sli == max(sli, na.rm=TRUE), arr.ind=TRUE)
      sli.max <- sli[k]
      thresh <- mean(c(bg, sli.max))
    } else {
      thresh <- mean(c(bg, suv.max))
    }
    vol.mask[, , z[i]] <- sli >= thresh & mask[, , z[i]]
    vol.mask[is.na(vol.mask)] <- FALSE
    vol[i] <- sum(vol.mask[, , z[i]], na.rm=TRUE)
    suv.total[[i]] <- sli[vol.mask[, , z[i]]]
  }
  list(totalVolume=sum(vol), totalSUV=mean(unlist(suv.total)),
       mask=vol.mask)
}

