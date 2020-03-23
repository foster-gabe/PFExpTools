#' fullCurate
#'
#' Are you feeling really, spectacularly lazy? This is the function for you. If
#' you know the GEO accession number of your expression set and the columns
#' necessary for nameChange, this function will do the entire curation for you
#' in one call.
#'
#' @param GSE This is the GSE accession number.
#' @param platcols These are the columns necessary for the nameChange function;
#'   see nameChange for details.
#' @param method These are the available methods for name changes; see
#'   nameChange for details.
#' @param aliases IF method is 'aliases', list the alias file name or 'current'
#'   for a PlasmoDB download. See nameChange for details.
#' @param transcripts If method is 'transcripts', this is the transcriptome
#'   fasta file, or if 'current' it will download the most recent file.
#' @param match see curateExpData.
#' @param secmatch See curateExpData.
#' @param pct See curateExpData.
#' @param pullmeta See curateExpData.
#'
#' @return See curateExpData.
#' @export
#'
#' @examples
fullCurate <- function(GSE, platcols = c(1,2), method = "blast", aliases = "current",
                       transcripts = "current", match = 130, secmatch = 60,
                       pct = 0.8, pullmeta = T){

  alldata <- GEOquery::getGEO(GEO = GSE, destdir = getwd())

  platform <- Biobase::annotation(alldata[[1]])

  namelist <- nameChange(platform = platform,
                         platcols = platcols,
                         method = method,
                         aliases = aliases,
                         transcripts = transcripts)

  output <- curateExpdata(GSE,namelist,local = F,pct = pct,pullmeta = pullmeta)

  output


}
