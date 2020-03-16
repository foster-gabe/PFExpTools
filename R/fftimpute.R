#' fftimpute
#'
#' This function imputes timepoints using the Fourier Transform; strongly recommend using another option here, as this is INCREDIBLY specific to a 2 hr timecourse over 56 hr
#' @param expdata An expression matrix with genes as rows and samples as columns
#' @param times A matrix with row names as samples and one column containing sampling times
#'
#' @return This function returns an expression matrix containing both the measured and imputed timepoints
#' @export
#'
#' @examples
fftimpute <- function(expdata, times){


  firstsample <- min(times)

  fulltimes <- firstsample:(max(times)+firstsample-1)


  expdata <- rbind(t(times), expdata)
  rownames(expdata)[1] <- "Time"

  #build list of new timepoints to impute
  newtimes <- subset(fulltimes,  !(fulltimes %in% times))

  #build matrix of new timepoints to impute
  spacer <- matrix(data = NA, ncol = length(newtimes), nrow = nrow(expdata))

  rownames(spacer) <- rownames(expdata)



  spacer[1,] <- newtimes



  #bind new times to existing data, sort by Time
  fulltimecourse <- cbind(expdata, spacer)
  fulltimecourse <- fulltimecourse[,order(fulltimecourse[1,])]
  fulltimecourse <- apply(fulltimecourse, 1, unlist)

  for(i in 1:ncol(fulltimecourse)){


    origfft <- fft(unname(expdata[i,]))

    numtimepoints <- ncol(expdata)

    if(numtimepoints %% 2 == 0){

      spacedfft <- c(origfft[1:(numtimepoints/2)], rep(0,numtimepoints), origfft[((numtimepoints/2)+1):numtimepoints])
    }else{

      spacedfft <- c(origfft[1:(floor(numtimepoints)/2)], rep(0,numtimepoints), origfft[ceiling(numtimepoints/2):numtimepoints])

    }


    imputedcourse <- fft(spacedfft, inverse = T)
    imputedcourse <- imputedcourse / length(expdata[i,])
    fulltimecourse[,i] <- Re(imputedcourse)


  }


  nameframe <- rownames(fulltimecourse)
  nameframe <- as.data.frame(nameframe, stringsAsFactors = F)
  nameframe[,2] <- fulltimes
  prefix <- stringr::str_extract(nameframe[1,1], "[A-Z]*[0-9]*[A-Z]*_")
  suffix <- stringr::str_extract(nameframe[1,1], "[0-9]*_1")


  for(i in 1:nrow(nameframe)){
    nameframe[i,1] <- paste(prefix, nameframe[i,2], "hpi", suffix, sep = "")
  }
  rownames(fulltimecourse) <- nameframe$nameframe

  fulltimecourse <- t(fulltimecourse)

  fulltimecourse

}
