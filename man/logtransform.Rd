% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LogTransform.R
\name{LogTransform}
\alias{LogTransform}
\title{LogTransform}
\usage{
LogTransform(newdata, refavg)
}
\arguments{
\item{newdata}{This is the data to be normalized; this must be a numeric matrix with rows as genes and columns as samples.}

\item{refavg}{This is the list of gene averages across a time course; must be a data frame with 1 column, with row names as genes.}
}
\value{
This function returns a matrix in an identical format to newdata, with all cells being log2 ratio transformed.
}
\description{
If your microarray data isn't log2 ratio transformed, use this. This will takeyour expression data and log2 ratio normalize it to a set of gene averages you provide.
}
