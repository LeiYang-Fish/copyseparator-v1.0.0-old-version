#' @title copy_validate
#'
#' @description This function checks if chimeric sequences may have been formed during the assembling of gene copies.
#'
#' @param filename A DNA alignment in fasta format that contains sequences of two or more gene copies.
#'
#' @param copy_number An interger (e.g. 2,3, or 4) giving the number of gene copies in the input file.
#'
#' @param read_length An interger (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data.
#'
#' @return A histogram in pdf format showing the comparison between the distance among neighboring variable sites and the read length.
#'
#' @importFrom seqinr read.fasta
#'
#' @importFrom grDevices dev.off pdf
#'
#' @importFrom graphics abline hist text
#'
#' @importFrom beepr beep
#'
#' @export copy_validate
#'

copy_validate<-function(filename,copy_number,read_length)
{
  sink("log.txt", append=FALSE, split=TRUE) # begin to record log
  error_log_function <- function() {
    cat(geterrmessage(), file="Error_log.txt", append=T) # begin to record error log
  }

  all_copies <- seqinr::read.fasta(filename, seqtype = "DNA", as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE)

  all_dist <-  character(0)
  for (i in 1:(copy_number)) {
    for (j in i:(copy_number)) {
      Paralog_diff <- which(unlist(strsplit(all_copies[[i]], ""))!=unlist(strsplit(all_copies[[j]], "")))

      dist_var_sites <- sapply(2:length(Paralog_diff), function(x) Paralog_diff[x]-Paralog_diff[x-1])
      all_dist <- append(all_dist,  dist_var_sites)
      cat(paste0("Copy ", i, " and Copy ", j, " is different at site: ", Paralog_diff,"\n"))
    }
  }

  dist_final <- sort(as.numeric(all_dist))

  if (max(dist_final)>read_length)
   cat ("Warning! Distance between variable sites is larger than the read length (see the plot in pdf), chimeric sequences may have formed during assembling!", "\n")

  text_to_show <- "read_length"

  grDevices::pdf(file="Variable sites distribution among gene copies.pdf", width=8, height=8)

  graphics::hist(dist_final,main="Distance among neighboring variable sites VS. Read length", xlab="Base pairs",
       ylim=c(0,length(dist_final)/2),xlim=c(0,read_length+100),breaks=20,col="blue")

  graphics::abline(v=read_length, col="red", lwd=3, lty=2)
  graphics::text(read_length, length(dist_final)/2, text_to_show, pos = 4, offset = 0.5)

  beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
  options("error" = error_log_function)
  sink() # turn off log
  grDevices::dev.off()
}




