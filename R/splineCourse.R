#' splineCourse
#'
#' This function takes an existing time course and uses spline interpolation to
#' create an imputed time course with the desired resolution.
#'
#' @param expdata Actual expression data; a matrix with genes as rows and time
#'   course samples as columns. In order from earliest to latest.
#' @param times The actual sample times for the expression matrix.
#' @param fulltimes A vector containing the times desired; for example, an hourly
#'   time course over a single IDC would be 1:48.
#'
#' @return Returns an expression matrix expanded with additional columns to contain samples at the desired times.
#' @export
#'
#' @examples
splineCourse <- function(expdata, times, fulltimes = 1:48){

  output <- matrix(nrow = nrow(expdata), ncol = length(fulltimes))
  rownames(output) <- rownames(expdata)
  colnames(output) <- paste("hpi_", fulltimes, "")

  for(i in 1:nrow(expdata)){
    splinefit <- spline(y=expdata[i,],x=times, xout = fulltimes)
    output[i,] <- splinefit$y
  }
  output
}
