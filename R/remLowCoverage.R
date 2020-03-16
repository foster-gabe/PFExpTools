#' remLowCoverage
#'
#' Removes any gene from the expression dataset that does not have values in a minimum percentage of samples.
#'
#' @param expdata This is a matrix containing expression values' rows must be genes and columns must be samples.
#' @param mincov This is the minimum percent coverage; probes must appear in at least this percentage of samples to pass QC.
#' @param faillist Default is FALSE; if TRUE, function returns a list containing both the curated expression matrix and a list of genes that failed QC.
#'
#' @return This function returns an expression matrix with the genes that do not pass QC removed.
#' @export
#'
#' @examples
remLowCoverage <- function(expdata, mincov, faillist = F){
  failedprobes <- as.list("spacer")

  #generate list of probes not appearing in enough samples
  for(i in 1:nrow(expdata)){
    if(sum(is.na(expdata[i,]))/(ncol(expdata)) > 1-mincov){
      failedprobes <- c(failedprobes, rownames(expdata)[i])
    }

  }

  #remove failed probes
  expdata <- expdata[!rownames(expdata) %in% failedprobes,]
  failedprobes <- failedprobes[-1]
  if(faillist){return(list(expdata, failedprobes))}else{return(expdata)}

}






