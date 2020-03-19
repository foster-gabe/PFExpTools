#' peakPhase
#'
#' @param expdata This must be an expression matrix, with rows as genes and
#'   columns as samples. Sample names (column names) must match sample names in
#'   the "times" parameter.
#' @param times This must be a matrix with one column; the column is actual
#'   sampling time, and the row names are the sample names.
#' @param method This allows the user to choose the method for peak time
#'   detemination; current options are "poly"(default), "sin", "spline", and
#'   "FFT". Note, FFT will only work on evenly sampled time courses (or evenly
#'   imputed time courses).
#'
#' @return This function returns a data frame with a single column; the row
#'   names are the gene names, and the column contains the time at which the
#'   gene is calculated to express maximally.
#' @export
#'
#' @examples
peakPhase <- function(expdata, times, method = "poly"){

  if(method == "poly"){
    output <- polyPhase(expdata, times)
  }else if(method == "sin"){
    output <- sinPhase(expdata, times)
  }else if (method == "FFT"){
    output <- FFTphase(expdata, times)
  }else if (method == "spline"){
    output <- splinePhase(expdata, times)
  }else {stop("Invalid method chosen, use only 'poly', 'sin', 'spline', or 'FFT'")}

  output

}

#' polyPhase
#'
#' This is a sub function that calculates peak expression by 6th order
#' polynomial fit. Accessed through peakPhase.
#'
#' @param expdata see peakPhase
#' @param times see peakPhase
#'
#' @return see peakPhase
#'
#' @examples
polyPhase <- function(expdata, times){

  fitfunc<-function(x,y){
    y[1] + x*y[2] + x^2*y[3] + x^3*y[4] + x^4*y[5] + x^5*y[6] + x^6*y[7]
  }

  output <- data.frame(Phase = double())
  for(i in 1:nrow(expdata)){
    gene <- expdata[i,]
    polyfit <- lm(gene ~ poly(times,6, raw = T))
    #find maximum from formula
    max <- as.numeric(optimize(fitfunc, 1:48, maximum = T, y = polyfit$coefficients))
    output[i,1] <- max[1]
    rownames(output)[i] <- rownames(expdata)[i]
  }
  output

}

#' sinPhase
#'
#' This is a sub function that calculates peak expression by sin fitting.
#' Accessed through peakPhase.
#'
#' @param expdata see peakPhase
#' @param times see peakPhase
#'
#' @return see peakPhase
#'
#' @examples
sinPhase <- function(expdata, times){

  output <- data.frame(Phase = double())
  for(i in 1:nrow(expdata)){

    gene <- expdata[i,]

    xc <- cos(2*pi*times/48)
    xs <- sin(2*pi*times/48)

    sinfit <- lm(gene ~ xc + xs, na.action = na.exclude)
    imputegene <- predict(sinfit)
    fitfunc<-function(x){sinfit$coefficients[2]*cos(2*pi*x/48) + sinfit$coefficients[3]*sin(2*pi*x/48) + sinfit$coefficients[1]}
    max <- as.numeric(optimize(fitfunc, 1:48, maximum = T)[1])
    output[i,1] <- max[1]
    rownames(output)[i] <- rownames(expdata)[i]

    }
  output

}

#' FFTphase
#'
#' This function calculates peak phase for each gene in a time course using the
#' Fourier Transform.
#' @param expdata see peakPhase
#' @param times see peakPhase
#'
#' @return see peakPhase
#'
#' @examples
FFTphase <- function(expdata, times){
  output <- data.frame(Phase = double())
  for(i in 1:nrow(expdata)){
    genefft <- fft(expdata[i,])
    maxphase <- atan2(-Im(genefft[2]), Re(genefft[2]))*(24/pi)
    output[i,1] <- maxphase
    if(output[i,1] < 0){output[i,1] <- output[i,1] + 48}
    rownames(output)[i] <- rownames(expdata)[i]
  }
  output

}

#' splinePhase
#'
#' This function fits a spline to the expression data and determines the time at
#' which each gene is maximally expressed. Called by peakPhase.
#'
#' @param expdata see peakPhase
#' @param times see peakPhase
#'
#' @return see peakPhase
#'
#' @examples
splinePhase <- function(expdata, times){

  output <- data.frame(Phase = double())
  for(i in 1:nrow(expdata)){
    gene <- expdata[i,]
    splinefit <- smooth.spline(times, gene, df=15)
    imputegene <- predict(splinefit, x=seq(0,48,.1))
    max <- imputegene[[1]][which.max(imputegene[[2]])]
    output[i,1] <- max[1]
    rownames(output)[i] <- rownames(expdata)[i]

  }

  output

}
