#' shuffle reads
#' 
#' Uses the uShuffle library to shuffle reads
#' 
#' @param reads An object of \link[Biostrings:XStringSet-class]{BStringSet}.
#' @param k the k-let size.
#' @param n the number of random sequences to generate.
#' @return An object of \link[Biostrings:XStringSet-class]{BStringSet}.
#' @import Biostrings
#' @importFrom methods is
#' @useDynLib uShuffleR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export 
#' @examples 
#' f <- system.file("extdata", "test.fa", package="uShuffleR")
#' shuffle(f)

shuffle <- function(reads, k=2, n=2){
  if(is.character(reads)){
    reads <- readDNAStringSet(reads)
  }
  stopifnot("reads must be an object of BStringSet"=
              inherits(reads, c("BStringSet", "DNAStringSet", "RNAStringSet",
                                "AAStringSet", "DNAString", "RNAString", 
                                "AAString")))
  if(is(reads, "DNAString")){
    reads <- DNAStringSet(reads)
  }
  if(is(reads, "RNAString")){
    reads <- RNAStringSet(reads)
  }
  if(is(reads, "AAString")){
    reads <- AAStringSet(reads)
  }
  if(length(reads)==0){
    return(NULL)
  }
  in_seqs <- as.character(reads)
  out_seqs <- rushuffle(in_seqs, k, n)
  if(is(reads, "DNAStringSet")){
    seqs <- DNAStringSet(out_seqs)
  }
  if(is(reads, "RNAStringSet")){
    seqs <- RNAStringSet(out_seqs)
  }
  if(is(reads, "AAStringSet")){
    seqs <- AAStringSet(out_seqs)
  }
  if(length(names(reads))){
    names(seqs) <- paste0(rep(names(reads), each=n), 
                          "_shuffle_", 
                          rep(seq.int(n), length(reads)))
  }
  return(seqs)
}