#' CRCgeneration
#'
#' @param stages This is a single column matrix; row names are sample names, and
#'   the column contains the calculated IDC location for each sample (in hpi)
#'
#' @return This function returns a 3 column matrix; the first column contains
#'   the provided location in hpi, and the additional columns contain the
#'   calculated x and y covariates, respectively.
#' @export
#'
#' @examples
CRCgeneration <- function(stages){

  stages <- as.data.frame(stages)
  stages[,2] <- cos((2*pi*stages[,1])/48)
  stages[,3] <- sin((2*pi*stages[,1])/48)

  colnames(stages) <- c("stagedTime", "xcov", "ycov")
  stages <- as.matrix(stages)
  stages

}
