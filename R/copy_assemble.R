#' @title copy_assemble
#'
#' @description Assembles a small number of overlapping DNA sequences into their respective gene copies.
#'
#' @param filename A fasta alignment of a small number of overlapping DNA sequences (results from "copy_separate") covering the entire length of the target gene. Check the alignment carefully before proceeding.
#'
#' @param copy_number An integer (e.g. 2,3, or 4) giving the anticipated number of gene copies. Must be the same value as used for "copy_separate".
#'
#' @param verbose Turn on (verbose=1; default) or turn off (verbose=0) the output.
#'
#' @return A fasta alignment of the anticipated number of full-length gene copies.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_order
#'
#' @importFrom Biostrings readDNAStringSet
#'
#' @importFrom DECIPHER ConsensusSequence
#'
#' @importFrom beepr beep
#' 
#' @examples 
#' \dontrun{
#' copy_assemble("inst/extdata/combined_con.fasta",2,1)
#' }
#'
#' @export copy_assemble
#'

copy_assemble<-function(filename,copy_number, verbose=1)
{
  sink("log.txt", append=FALSE, split=TRUE) # begin to record log
  error_log_function <- function() {
    cat(geterrmessage(), file="Error_log.txt", append=T)
  }

  if (copy_number<=1) stop ("The anticipated copy number must be a number larger than one!")

  Consensus_seq <- seqinr::read.fasta(file = filename,seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

  subset_sum <- length(Consensus_seq)/copy_number

  seq_to_con <- character(0)
  for (l in 1:copy_number) {
    for (i in seq(l, (subset_sum-1)*copy_number, copy_number)){
      ambiguity_num <- character(0)
      for (j in 1:copy_number) {
        seqinr::write.fasta(sequences = c(Consensus_seq[i], Consensus_seq[i+copy_number+j-l]), names(c(Consensus_seq[i], Consensus_seq[i+copy_number])),file.out = paste0("cnx_",i,"_",i+copy_number+j-l,".fasta"))
        ambiguity_num <- append(ambiguity_num, sum(as.integer(stringr::str_count(as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0("cnx_",i,"_",i+copy_number+j-l,".fasta"), format="fasta",nrec=-1L, skip=0L),
                                                                                                       threshold = 0.4,ambiguity = TRUE, noConsensusChar = "N")), c("M", "K", "R", "Y","N","W", "S", "H", "V", "D", "B")))))
        if (verbose) { cat(paste0("Matching sequences ",i, " & ", i+copy_number+j-l, "\n"))}
      }
      if (verbose) { cat("Number of ambiguous sites for each match:", ambiguity_num, "\n")}
      if (verbose) { cat(paste0("--- The pair ", which(as.numeric(ambiguity_num)==min(as.numeric(ambiguity_num))), " has fewer ambiguous sites and should be assembled together\n"))}

      seq_to_con <- append(seq_to_con, paste0("A",i,"B_A",i+copy_number+which(as.numeric(ambiguity_num)==min(as.numeric(ambiguity_num)))-l,"B"))
    }
  }
  seq_to_con1 <- seq_to_con[stringr::str_order(seq_to_con, decreasing = FALSE, na_last = TRUE, locale = "en",numeric = TRUE)] # order numerically

  unlink("cnx_*") # Delete all intermediate files whose names begin with "cnx_"

  if (verbose) { cat("*************************************************************************\n")}
  # To find out which sequences belong to which gene copy
  Copy_list <- data.frame()
  for (i in 1:copy_number) {
    cp <- c(i)
    for (j in 1:(subset_sum-copy_number)) {
      m <- which(gsub(".*_", "", seq_to_con1[cp[j]]) == gsub("_.*", "", seq_to_con1)) # From A1B_A4B to A4B_A5B, then to A5B_A8B, then to A8B_A10B ...; linking 1, 4, 5, 8, 10 ...
      cp <- append(cp,m)
    }
    cp_all <- as.numeric(gsub("[AB]","",c(gsub("_.*", "", seq_to_con1[cp]), gsub(".*_", "", seq_to_con1[cp][subset_sum-copy_number+1]))))
    Copy_list <- append(Copy_list,as.data.frame(cp_all))
  }

  seq_list <- character(0)
  for (i in 1:copy_number) {
    cat("List of sequences to be assembled for gene copy", i, ": ", Copy_list[[i]], "\n")
    for (j in (1:copy_number)[-i]){
      seq_list <- append(seq_list, intersect(as.character(Copy_list[[i]]),as.character(Copy_list[[j]])))
    }
  }

  cat("--- Sequences involved in the assembling of multiple gene copies: ", stringr::str_sort(unique(seq_list), numeric=TRUE), "\n")
  cat("Warning! If there are sequences involved in the assembling of multiple gene copies, please check your input file carefully and try to do the assembling again or do it manually!\n")


  all_copies_final <- character(0)
  for (i in 1:copy_number) {
    copy_subseqs <- lapply(as.numeric(unlist(Copy_list[i])), function(x) Consensus_seq[x])
    copy_subseqs_name <- lapply(as.numeric(unlist(Copy_list[i])), function(x) names(Consensus_seq[x]))
    seqinr::write.fasta(sequences = copy_subseqs, copy_subseqs_name,file.out=paste0("Copy_",i,"_subseqs.fasta"))
    seqinr::write.fasta(sequences = as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0("Copy_",i,"_subseqs.fasta"), format="fasta",nrec=-1L, skip=0L),
                                                           threshold = 0.4,ambiguity = TRUE, noConsensusChar = "N")), paste0("Copy_",i,"_final"),file.out=paste0("Copy_",i,"_final.fasta"))
    all_copies_final <- append(all_copies_final, seqinr::read.fasta (paste0("Copy_",i,"_final.fasta"), seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE))
  }
  seqinr::write.fasta(sequences = all_copies_final, names(all_copies_final),file.out=paste0("All_final_copies.fasta"))
  unlink("*_final.fasta")

  if (verbose) { cat("*************************************************************************\n")}
  cat("Run finished!\n")
  beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
  options("error" = error_log_function)
  sink() # turn off log
}
