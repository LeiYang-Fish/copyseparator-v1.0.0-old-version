#' toydata
#' 
#' Provides small example dataset for major functions. They are located at: /inst/extdata
#'
#' The file "toydata.fasta" is for functions "copy_separate" and "subset_downsize".
#' The file "combined_con.txt" is for the function "copy_assemble".
#' The file "toysubset.fasta" is for the function "copy_detect".
#' The file "Final_two_copiess.fasta" is for the function "copy_validate".
#'
#' @docType data
#' @usage
   # copy_separate("inst/extdata/toydata.fasta",2,300,225,1)
   # copy_assemble("inst/extdata/combined_con.fasta",2,1)
   # copy_validate("inst/extdata/Final_two_copies.fasta",2,300,1)
   # subset_downsize("inst/extdata/toydata.fasta",300,225,1)
   # copy_detect("inst/extdata/toysubset.fasta",2,1)

### After the test runs, you may run the following codes to delete all the resulting files.
   # unlink("*.fasta")
   # unlink("*.txt")
   # unlink("*.pdf")
   # unlink("inst/extdata/toydata_intermediate_files", recursive = TRUE)
   # unlink("inst/extdata/toydata_combined*.txt")
   # unlink("inst/extdata/toysubset_*")

#' @format fasta
#' @references This data set was collected by the author and is unpublished.
#' @keywords dataset
