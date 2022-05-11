#' toydata
#'
#' Provides toydata (in the folder: "data") only for package testing during package submission. The real example dataset is located at: /inst/extdata.
#'
#' The file "toydata.fasta" is for functions "copy_separate" and "subset_downsize".
#' The file "combined_con.txt" is for the function "copy_assemble".
#' The file "toysubset.fasta" is for the function "copy_detect".
#' The file "Final_two_copiess.fasta" is for the function "copy_validate".
#'
#' @docType data
#' @format fasta
#' @references This data set was collected by the author and is unpublished.
#' @keywords dataset
#' @examples
#'
input_data1 <- get(load("data/toydata.RData"))
seqinr::write.fasta(sequences = input_data1, names = names(input_data1), file.out = "data/toydata.fasta")
input_data2 <- get(load("data/combined_con.RData"))
seqinr::write.fasta(sequences = input_data2, names = names(input_data2), file.out = "data/combined_con.fasta")
input_data3 <- get(load("data/Final_two_copies.RData"))
seqinr::write.fasta(sequences = input_data3, names = names(input_data3), file.out = "data/Final_two_copies.fasta")
input_data4 <- get(load("data/toysubset.RData"))
seqinr::write.fasta(sequences = input_data4, names = names(input_data4), file.out = "data/toysubset.fasta")


#' copy_separate
copy_separate("data/toydata.fasta",2,300,225,1)

#' copy_assemble
copy_assemble("data/combined_con.fasta",2,1)

#' copy_validate
copy_validate("data/Final_two_copies.fasta",2,300,1)

#' subset_downsize
subset_downsize("data/toydata.fasta", 300,225,1)

#' copy_detect
copy_detect("data/toysubset.fasta",2,1)



#' delete all resulting files after testing
unlink("*.fasta")
unlink("*.txt")
unlink("*.pdf")
unlink("data/toydata_intermediate_files", recursive = TRUE)
unlink("data/toydata_combined*.txt")
unlink("data/toysubset_*")
unlink("data/*.fasta")




