#' LogTransform
#'
#' If your microarray data isn't log2 ratio transformed, use this. This will takeyour expression data and log2 ratio normalize it to a set of gene averages you provide.
#' @param newdata This is the data to be normalized; this must be a numeric matrix with rows as genes and columns as samples.
#' @param refavg This is the list of gene averages across a time course; must be a data frame with 1 column, with row names as genes.
#'
#' @return This function returns a matrix in an identical format to newdata, with all cells being log2 ratio transformed.
#' @export
#'
#' @examples
LogTransform <- function(newdata,refavg){
  #sort data and reference averages

  if(!isTRUE(all.equal(rownames(newdata), rownames(refavg)))){stop("Gene names do not match")}

  newdata <- newdata[order(rownames(newdata)),]
  refavg <- as.matrix(refavg[order(rownames(refavg)),])

  for (i in 1:nrow(newdata)){
    for (j in 1:ncol(newdata)){

      newdata[i,j] <- log2((abs(newdata[i,j])^sign(newdata[i,j])) / refavg[i])
    }
  }

  newdata

}



#' GetAverages
#' This function will get the average expression of all genes across a reference time course. Note: the reference time course should be evenly samples across a single IDC.
#' @param reftc This is a reference time course to average; this must be a numeric matrix with genes as rows and samples as columns.
#'
#' @return This function returns a matrix with row names as gene names and a single column containing the average expression across the reference time course.
#' @export
#'
#' @examples
getAverages <- function(reftc){
  output <- as.matrix(apply(reftc, 1, mean, na.rm = T))
  output
}


