#' Provides example dataset for major functions. They are located at: /inst/extdata
#'
#' The example file "Sucker.fasta" is for functions "copy_separate" and "subset_downsize". It is the results of 69,508 Illumina reads mapped to a reference gene (length: ~1000bp).
#' The data were obtained by adding adapters to the PCR amplicons and sent for sequencing on the illumina Miseq platform (2x300).
#' The reference sequence and reads mapped outside the boundaries defined by both PCR primers have been removed after the mapping using Geneious.
#' The example file "Sucker_combined_consensus_2copies_overlap225.txt" is for the function "copy_assemble". It is the resulting file from "copy_separate".
#' The example file "Subset_1_downsized.fasta" is for the function "copy_detect". It is one of the 12 subsets from the results of "copy_separate".
#' The example file "All_final_copies.fasta" is for the function "copy_validate". It contains the two gene copy sequences assembled by "copy_assemble".
#'
#' @docType data
#'
#' @usage
#' copy_separate("Sucker.fasta",2,300,225)
#'
#' subset_downsize("Sucker.fasta", 300,225)
#'
#' copy_detect("Subset_1_downsized.fasta",2,0.5)
#'
#' copy_assemble("Sucker_combined_consensus_2copies_overlap225.txt",2)
#'
#' copy_validate("All_final_copies.fasta",2,300)
#'
#' @format fasta
#' @references This data set was collected by the author and is unpublished.
#' @keywords dataset
#'
#'
#'
#'
