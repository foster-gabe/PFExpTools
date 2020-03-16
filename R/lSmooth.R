#' lSmooth
#'
#'This function will loess smooth a time course, first demonstrated on P. fal time courses by Bozdech et al 2003
#' @param expdata Expression data matrix, with columns as samples and rows as genes
#' @param time matrix with one column; names are samples and column contains sample time
#' @param smoothvalue loess smoothing value; default is 0.3
#'
#' @return
#' @export
#'
#' @examples
lSmooth <- function(expdata, time, smoothvalue = 0.3){
  if(!isTRUE(all.equal(colnames(expdata), rownames(time)))){stop("Sample names do not match")}

  output <- expdata
  expdata <- as.data.frame(expdata)

  for(i in 1:nrow(expdata)){
    test <- loess(unlist(expdata[i,]) ~ time, data = expdata, span = smoothvalue, na.action = na.exclude)
    smoothed <- predict(test)
    output[i,] <- smoothed
  }
  output
}





