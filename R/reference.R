##
## Copyright (c) 2012, Brandon Whitcher
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

simplifiedReferenceTissueModel <- function(conc, time, mask, k2ref=NULL) {
  require("minpack.lm")
  if (is.null(k2ref)) {
    func.model <- 
  } else {
    
  }
  func <- function(theta, signal, time, ...) {
    out <- signal - func.model(time, theta, p)
    out[!is.na(out)]
  }
  nvoxels <- sum(mask)
         weinmann = ,
         weinmann.empirical = ,
         kety.orton.exp = ,
         kety.orton.cos = {
           if (is.null(guess)) {
             guess <- c("th1"=-1, "th3"=-1)
           } else {
             if (length(guess) != 2 || ! all(names(guess) %in% c("th1","th3"))) {
               stop("Names of starting parameters must be \"th1\" and \"th3\"")
             }
           }
         },
       )
  return(1)
}

multilinearReferenceTissueModel <- function(nim, mask, k2ref=NULL) {
  return(1)
}
