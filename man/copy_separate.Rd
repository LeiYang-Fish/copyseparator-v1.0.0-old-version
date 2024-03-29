% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copy_separate.R
\name{copy_separate}
\alias{copy_separate}
\title{copy_separate}
\usage{
copy_separate(filename, copy_number, read_length, overlap = 225, verbose = 1)
}
\arguments{
\item{filename}{A fasta file contains thousands of short reads that have been mapped to a reference. The reference and reads that are not directly mapped to the reference need to be removed after mapping.}

\item{copy_number}{An integer (e.g. 2,3, or 4) giving the anticipated number of gene copies in the input file.}

\item{read_length}{An integer (e.g. 250, or 300) giving the read length of your Next-generation Sequencing data. This method is designed for read length >=250bp.}

\item{overlap}{An integer describing number of base pairs of overlap between adjacent subsets. More overlap means more subsets. Default 225.}

\item{verbose}{Turn on (verbose=1; default) or turn off (verbose=0) the output.}
}
\value{
A fasta alignment of a small number of overlapping DNA sequences covering the entire length of the target gene. Gene copies can be assembled by reordering the alignment manually or use the function "copy_assemble".
}
\description{
Separates two or more gene copies from short-read Next-Generation Sequencing data into a small number of overlapping DNA sequences.
}
\examples{
\dontrun{
copy_separate("inst/extdata/toydata.fasta",2,300,225,1)
}

}
