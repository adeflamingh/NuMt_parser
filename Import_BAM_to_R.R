library(Rsamtools)
library(GenomeInfoDb)

#change to working directory
setwd("~/Box Sync/aDNA/lions (deflami2@illinois.edu)/NUMTS/Analysis/R_analysis")


##LOAD MITO AND NUMT BAM FILES##


#function for reading bamfile - from <script src="https://gist.github.com/SamBuckberry/9914246.js"></script>

readBAM <- function(bamFile){
  
  bam <- scanBam(bamFile)
  
  # A function for collapsing the list of lists into a single list
  # as per the Rsamtools vignette
  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  #return a list that can be called as a data frame
  return(bam_df)
}

# specify and load the bam file
bamFile <- "Singer1_filtered.bam" #path to mitochondrial bam file
mitobam <- readBAM(bamFile)

bamFile <- "Singer1_NuMt.bam" #path to numt bam file
numtbam <- readBAM(bamFile)

