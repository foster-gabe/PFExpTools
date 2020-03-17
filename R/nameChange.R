
#' nameChange
#'
#' This function takes a platform file from GEO (a cut and paste from the full
#' table view, really), and identifies the new genome assembly name for each
#' probe using one of a few methods. This function is a wrapper for subfunctions
#' that carry out each method.
#'
#' @param platform This is the filename for the platform from GEO. A simple
#'   cut/paste in to a txt file will work here (it's tab delimited by default).
#' @param method There are 3 methods here; "swap" just takes the new gene name
#'   from the platform file (least accurate), "alias" will find new names based
#'   on the current PlasmoDB assembly file (not bad), "blast" will take a
#'   provided Annotated Transcripts file from PlasmoDB, blast each probe against
#'   it, and only keep probes that meet the bitscore criteria (default is 130
#'   for perfect hit and 60 for minimal secondary hits).
#' @param aliases This is the flat alias file created by the function
#'   "flattenAlias"; only necessary if selected method is "alias".
#' @param transcripts The fasta file containing the annotated transcripts; this
#'   can be obtained from PlasmoDB.
#' @param match The bitscore for a perfect primary match; only necessary if
#'   "blast" is selected, default is 130.
#' @param secmatch The maxmium bitscore for secondary probe alignments; only
#'   necessary if "blast is selected, default is 60.
#'
#' @return
#' @export
#'
#' @examples
nameChange <- function(platform, method = "swap", aliases = NA, transcripts = NA, match = 130, secmatch = 60){
  #read and curate platform
  suppressWarnings(platdata <- read.table(platform, skip = which.max(count.fields(platform)) - 1, sep = "\t",header = T, stringsAsFactors = F))
  #Remove blank row(s)
  platdata <- platdata[which(grepl("[a-zA-Z0-9]", platdata[,1])),]

  #Farm out to sub functions by method, catch obvious errors
  if(method == "swap"){
    output <- swapChange(platdata)
  } else if(method == "alias"){
    if(!all(is.na(aliases))){
      output <- aliasChange(platdata, aliases)
    }else{stop("No alias file provided for method alias")}
  }else if(method == "blast"){
    if(!is.na(transcripts)){
      output <- blastChange(platdata, platform, transcripts, match, secmatch)
    }else{stop("No transcript file provided for method blast")}
  }else{stop("Invalid method provided, only supported methods are swap, alias and blast")}

  output

}

#' aliasChange
#'
#' This function performs the alias based gene renaming; called by nameChange
#'
#' @param platdata see nameChange
#' @param aliases see nameChange
#'
#' @return see nameChange
#'
#' @examples
aliasChange <- function(platdata, aliases){

  #pull probe ID and old gene name
  output <- platdata[,c("ID", "ORF_old")]

  #remove probes with no gene name
  output <- output[which(nchar(output[,2]) > 1),]

  #create named vector from aliases for easy lookup
  aliasvector <- aliases[,2]
  names(aliasvector) <- aliases[,1]
  #Perform lookup
  output[,3] <- sapply(output[,2], function(x) {unname(aliasvector[x])})
  #remove any probes not mapping to new gene ID
  output <- output[which(!is.na(output[,3])),]
  #pull only probeID and new gene ID, name and return
  output <- output[,c(1,3)]
  colnames(output) <- c("Probe", "GeneID")
  output
}

#' swapChange
#'
#' This function simply pulls the new ORF and probe ID and creates a renaming
#' file from that; I swear this wasn't there before, but now that we have the
#' new GeneIDs listed might as well build a function to rename based on them.
#' This function is called by nameChange.
#'
#' @param platdata see nameChange
#'
#' @return see nameChange
#'
#' @examples
swapChange <- function(platdata){

  output <- platdata[,c("ID", "ORF")]
  colnames(output) <- c("Probe", "GeneID")
  output <- output[which(grepl("[a-zA-Z0-9]", output[,2])),]
  output
}

#' blastChange
#'
#' This function is...elaborate. It requires blast be installed and be in your
#' PATH. This function takes your platform data file and your annotated
#' transcriptome file and actually creates a blast database, searches all probes
#' against it, and identifies hits that are unique and maximal (depending on
#' your match and secondary match parameters.) It's called from nameChange.
#'
#' @param platdata see nameChange
#' @param platform curated data from nameChange
#' @param transcripts see nameChange
#' @param match see nameChange
#' @param secmatch see nameChange
#'
#' @return see nameChange
#'
#' @examples
blastChange <- function(platdata, platform, transcripts, match, secmatch){
  #make blast database
  system2(command = "makeblastdb", args = c("-in", transcripts, "-dbtype nucl"))

  #make fasta file from probeID/probe seqs
  probes <- platdata[,1:2]

  #first, add > to the beginning of probe lines
  probes[,1] <- gsub("^", ">", probes[,1])

  #check this out we can flip and flatten and it's perfect
  probes <- as.data.frame(t(probes), stringsAsFactors = F)
  probes <- as.data.frame(unlist(probes))
  filename <- paste(gsub(".txt", "", platform), ".fasta", sep = "")
  write.table(probes, file = filename, col.names = F, quote = F, row.names = F)

  #cool now we have a fasta file to search across, and we have a database built
  #so let's get searching shall we

  blastresults <- system2(command = "blastn", args = c("-query", filename,
                                                       "-db", transcripts,
                                                       "-outfmt", 10),
                                                       wait = T,
                                                       stdout = T)

  blastresults <- read.table(textConnection(blastresults), sep = ",")

  #and perfectly formatted blast results. Let's pull all probes that meet our
  #minimum match score
  blastresults <- blastresults[which(blastresults$V12 >= secmatch),]

  #now let's remove all probes that have multiple quality hits
  blastresults <- blastresults[which(!duplicated(blastresults$V1)),]

  #and now let's pull only perfect matches
  blastresults <- blastresults[which(blastresults$V12 >= match),]

  #we don't care about exons so let's remove the splicing indicators in the new
  #gene names
  blastresults[,2] <-sapply(blastresults[,2], function(x) {gsub("\\.\\d{1}$", "", x)})

  #homestretch, let's pull the probes and new IDs and return them

  output <- blastresults[,1:2]
  colnames(output) <- c("Probe", "GeneID")
  output
}


#' flattenAlias
#'
#' The alias file from PlasmoDB is unwieldy and awful. Here is a function for
#' flattening it for easier find and replace.
#'
#' @param aliasfilename File name of alias file obtained from PlasmoDB.
#'
#' @return Data frame with two columns; first column is old identifiers, and
#'   second column is their mapping to the most recent name (per the alias file)
#' @export
#'
#' @examples
flattenAlias <- function(aliasfilename){
  #file is \t delim, with uneven numbers of columns which is a joy. We'll create
  #a reverse flattened file here.
  aliases <- read.delim(aliasfilename, sep = " ", stringsAsFactors = F, header = F)
  output <- as.data.frame(matrix(nrow = 1, ncol = 2))
  for(i in 1:nrow(aliases)){
    line <- unlist(strsplit(aliases[i,1], "\t"))
    for(j in 2:length(line)){
      output <- rbind(output, c(line[j], line[1]))
    }
  }
  output <- output[which(!is.na(output[,1])),]

  #Our output is close; it's still kind of a hot mess with alternative splicing
  #and duplicates in here. Let's clean up as best we can.

  #First, we're not catching any alternate splicing with these old probes so
  #let's strip the splicing IDs
  output[,2] <-sapply(output[,2], function(x) {gsub("\\.\\d{1}$", "", x)})
  #collapse identical rows
  output <- dplyr::distinct(output)
  #some old genes are mapping to multiple new names; they need to be removed from our set.
  output <- output[which(!duplicated(output[,1])),]
  output
}



