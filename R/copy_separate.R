#' @title copy_separate
#'
#' @description This function separates two or more gene copies from short Next-Generation Sequencing reads.
#'
#' @param filename A fasta file containing thousands of short reads that have been mapped to a reference gene.
#'
#' @param copy_number An interger (e.g. 2,3, or 4) giving the expected number of gene copies in the input file.
#'
#' @param read_length An interger (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data. This method is designed for read length 250bp or longer.
#'
#' @param overlap An interger giving how many base pairs of overlap between adjacent subsets. Must be smaller than read_length. Default 225.
#'
#' @param lower_threshold A decimal (range: 0.01-0.99) giving the OTU identity cutoff for divisive heirarchical clustering analyses.Default 0.4.
#'
#' @param kmer An integer giving the k-mer size used to generate the input matrix for k-means clustering. Default 5.
#'
#' @return A fasta alignment of a small number of overlapping DNA sequences covering the entire length of the target gene. Gene copies can be assembled by moving sequences around.
#'
#' @importFrom seqinr read.fasta write.fasta
#'
#' @importFrom stringr str_count str_sort
#'
#' @importFrom kmer otu
#'
#' @importFrom Biostrings readDNAStringSet
#'
#' @importFrom ape as.character.DNAbin read.FASTA
#'
#' @importFrom DECIPHER ConsensusSequence
#'
#' @importFrom beepr beep
#'
#' @export copy_separate
#'

copy_separate<-function(filename,copy_number,read_length,overlap=225,lower_threshold=0.4, kmer=5)
{
sink("log.txt", append=FALSE, split=TRUE) # begin to record log
error_log_function <- function() {
  cat(geterrmessage(), file="Error_log.txt", append=T)
}

if (copy_number<=1) stop ("The expected copy number must be a number larger than one!")

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

# for each of the subset_1, 2, 3...
for (i in 1:length(Subsets)) {
  cat("********************************************\n")
  Subset <- ape::as.character.DNAbin(ape::read.FASTA(file=Subsets[i], type = "DNA"))
  cat((paste0("Clustering analyses for the Subset ", i,"\n")))

  # try different threshold values for OTU to find the major clusters (number=copy_number) for each subset
  for (j in seq(lower_threshold,0.99, by = 0.01)) {
    Subset_OTU <- kmer::otu(Subset, k = kmer, threshold = j, method = "central", nstart = 20)
    cat(paste0("threshold = ",j),"\n")
    cat(unique(Subset_OTU),"\n")
    if (length(unique(Subset_OTU))>=copy_number) {break}
  }
  reads_each_cluster <- sapply(unique(Subset_OTU), function(x) length(which(Subset_OTU==x)))

  cat(paste0("Best threshold = ",j),"\n")
  cat(unique(Subset_OTU),"\n")
  cat("Number of reads in each cluster\n")
  cat(reads_each_cluster,"\n")

  for (l in (1:copy_number)) {
    Picked_cluster <- Subset[which(Subset_OTU==unique(Subset_OTU)[which(reads_each_cluster==sort(reads_each_cluster)[length(unique(Subset_OTU))-l+1])])]
    seqinr::write.fasta(sequences = Picked_cluster, names = labels(Picked_cluster), file.out = paste0("Subset_",i,"_cluster_",l,".fasta"))
    cat(paste0("Number of reads in picked cluster ",l, " = ", length(Picked_cluster),"\n"))

    # calcuate the consensus sequence for the clusters of the subset
    seqinr::write.fasta(sequences = paste0(c(as.character(rep("-",(read_length-overlap)*(i-1))),as.character(DECIPHER::ConsensusSequence(Biostrings::readDNAStringSet(paste0("Subset_",i,"_cluster_",l,".fasta"),
                                                                                                                                        format="fasta",nrec=-1L, skip=0L),threshold = 0.4, ambiguity = TRUE, noConsensusChar = "N")[1])), collapse=''),
                names = paste0("Subset_",i,"_cluster_",l,"_consensus"), file.out = paste0("Subset_",i,"_cluster_",l,"_consensus.fasta"))
  }
}

# put together all the consensus sequences from all subsets into one file
Consensus_list <- stringr::str_sort(list.files(pattern="_consensus.fasta"), numeric = TRUE)
All_consensus <- lapply(1:length(Consensus_list), function (x) seqinr::read.fasta(file = Consensus_list[x], seqtype = "DNA",
                                                                          as.string = TRUE,forceDNAtolower = FALSE,set.attributes = FALSE, whole.header = TRUE))
filename_short <- gsub("[:.:].*","", filename) # remove file extensions, e.g. ".fasta", ".txt"

# move subset files into the intermediate files folder
dir.create(paste0(filename_short,"_intermediate_files"))
invisible(file.copy(list.files(pattern="Subset_"), paste0(filename_short,"_intermediate_files")))  # use "invisible" so that output do not show here
unlink(list.files(pattern="Subset_"))

seqinr::write.fasta(sequences=All_consensus, names=Consensus_list, file.out=paste0(filename_short,"_combined_consensus_",copy_number,"copies_overlap",overlap,".txt"))
cat("Run finished!\n")
beepr::beep(sound = 1, expr = NULL) # make a sound when run finishes
options("error" = error_log_function)
sink() # turn off log
}

