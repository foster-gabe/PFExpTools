#' tcImpute
#'
#' This function calls the bnstruct package "knn.impute" function to perform k nearest neighbor imputation on a time course to fill in missing data.
#'
#' @param expdata This is an expression matrix with rows as genes and columns as samples in a time course. Time points missing data should be inserted as columns in the correct place containing NAs.
#' @param neighbors This is the number of neighbors used in determining value.
#'
#' @return This fucntion returns the expression matrix with blank columns filled with imputed values.
#' @export
#'
#' @examples
tcImpute <- function(expdata, neighbors = 5){

  texpdata <- t(expdata)
  imputedata <- as.matrix(bnstruct::knn.impute(texpdata, k=neighbors))
  timputedata <- t(imputedata)
  timputedata


}
