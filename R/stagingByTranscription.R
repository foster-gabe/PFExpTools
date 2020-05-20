
#' stagingByTranscription
#'
#' This function compares each whole transcriptome sample to every time point in
#' the provided reference time course, and identifies the best time point fit.
#' Correlation is determined by the Pearson method.
#'
#' @param samples A matrix of expression samples, with rows as genes and columns
#'   as samples. This is the set that will be tested.
#' @param reference A reference expression time course, with rows as genes and
#'   columns as samples.
#' @param subset The default is to use all genes in common between the two sets;
#'   if a subset is desired, provide this as vector of gene names in the subset.
#' @param cortable If the full matrix of correlation values is desired, set to
#'   TRUE; default is FALSE.
#' @param heatmap A simple correlation heatmap will be generated as default. Set
#'   to false to skip.
#' @param start What hour the reference course starts at- default is 1hpi
#'
#' @return Default return is a matrix with one column; rows are samples and
#'   values are the time point of maximum correlation. The function will also
#'   plot a simple heatmap by default. If cortable is set to TRUE, the function
#'   will return a list; the first entry is the simple time point matrix, and
#'   the second entry is the full correlation matrix.
#' @export
#'
#' @examples
stagingByTranscription <- function(samples, reference, subset = NA, cortable = F, heatmap = T, start = 1){

  #Identify genes common in sample and reference set
  commongenes <- intersect(rownames(samples), rownames(reference))

  #If subset of staging genes is provided, use only those
  if(!is.na(subset)){commongenes <- intersect(commongenes, subset)}

  #reduce all sets to common genes
  samples <-  samples <- samples[match(commongenes,rownames(samples)),]
  reference <- reference <- reference[match(commongenes,rownames(reference)),]

  correlations <- cor(data.matrix(samples), data.matrix(reference), use="pairwise.complete.obs", method=c("pearson"))

  stages <- max.col(correlations) - (start - 1)
  stages <- as.data.frame(stages)
  row.names(stages) <- row.names(correlations)
  colnames(stages) <- "stagedTime"

  if(heatmap == T){
    heatmap(t(correlations), Colv=NA, Rowv=NA, labRow = T, xlab = "Sample Name", ylab = "Reference Timepoint (hpi)", margins = c(8,3))
  }

  if(cortable == T){

    output <- list("stages" = stages, "cortable" = correlations)

  } else {
    output <- stages
  }

  output

}
