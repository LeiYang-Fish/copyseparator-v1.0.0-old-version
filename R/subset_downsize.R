#' @title subset_downsize
#'
#' @description This function subdivides the imported read alignment into subsets and then downsizes each subset by deleting those sequences that have too many gaps or missing data.
#'
#' @param filename A fasta file containing thousands of short reads that have been mapped to a reference gene.
#'
#' @param read_length An interger (e.g. 250, or 300) giving the read length of your Next-generation sequencing data. This method is designed for read length 250bp or longer.
#'
#' @param overlap An interger giving how many base pairs of overlap between adjacent subsets. Must be smaller than "read_length".
#'
#' @return A number of overlapping subsets of the input alignment. More overlap means more subsets, but not fewer reads per subset.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_sort
#'
#' @importFrom beepr beep
#'
#' @export subset_downsize
#'
subset_downsize<-function(filename,read_length,overlap)
{
  sink("log.txt", append=FALSE, split=TRUE) # begin to record log
  error_log_function <- function() {
    cat(geterrmessage(), file="Error_log.txt", append=T)
  }

  if (read_length<250) warning ("This method is designed for read length 250bp or longer. Short reads can easily result in chimeric sequences.")

  NGS_reads <- seqinr::read.fasta(file = filename,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

  total_reads <- length(NGS_reads); cat(paste0("Total number of reads imported = ",total_reads,"\n"))
  alignment_length <- nchar(NGS_reads[[1]]); cat(paste0("Length of the alignment = ",alignment_length,"\n"))

  cat(paste0("Read length = ",read_length,"\n"))
  cat(paste0("Overlap between adjacent subsets = ",overlap,"\n"))
  if (overlap>=read_length) stop("Overlap between adjacent subsets must be smaller than the read length!")

  begin_number <- seq(1,alignment_length-200, by=read_length-overlap); cat("Beginning position of each subset","\n"); cat(begin_number,"\n")
  end_number <- begin_number+read_length-1; cat("Ending position of each subset\n"); cat(end_number,"\n")
  number_of_subsets <- length(begin_number); cat(paste0("Total number of subsets = ",number_of_subsets,"\n"))

  for (i in begin_number) {
    subset_original <- lapply(1:total_reads, function(x) {substr(NGS_reads[x],i,i+read_length-1)}) # begin to subdivide the big alignment into subsets, each has the length of read_length
    subset_small <- subset_original[which(as.character(lapply(1:total_reads, function(x) {substr(subset_original[x],1,10)}))!="----------")]
    subset_smaller <- subset_small[which(as.character(lapply(1:length(subset_small), function(x) {substr(subset_small[x],200,209)}))!="----------")]
    subset_smallest <- subset_smaller[which(lapply(1:total_reads, function(x) {stringr::str_count(substr(subset_smaller[x],1,200),"A")>25})==TRUE)]
    cat(paste0("Number of reads in subset ", which(begin_number==i), " = ",length(subset_smallest),"\n"))
    if (length(subset_smallest)<=2) stop("This subset has too few reads (<2). Enter a new value for the parameter 'overlap' and subset the dataset again!")
    seqinr::write.fasta(sequences = subset_smallest, names = 1:length(subset_smallest), file.out = paste0("Subset_",which(begin_number==i),"_downsized.fasta"))
  }

  Subsets <- stringr::str_sort(list.files(pattern="_downsized.fasta"), numeric = TRUE)
  cat("Subsetting and downsizing finished!\n")
  beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
}


