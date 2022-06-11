# How to install copyseparator
install.packages("copyseparator")
  
# Example on how to run copyseparator
https://github.com/LeiYang-Fish/copyseparator_more_examples

Due to the size limitation (25M) imposed by GitHub, only one example is provided in the above folder.

# What can copyseparator do?
Separate and assemble two or more gene copies from short-read Next-Generation Sequencing data for polyploids, gene families, mixed/contaminated samples etc.

# How many gene copies can copyseparator separate?
Two or more. However, it works better when there is fewer (e.g. two) gene copies. That is because copyseparator uses clustering analyses to separate reads from different gene copies. The more the gene copies, the harder it is to separate them apart.

# What is the length of gene copies that copyseparator can separate?
Theoretically, if the coverage is high and relatively even, there is no limitation on the length. In reality, however, copyseparator is usually used to separate gene copies that are several hundred to several thousand base pairs long. I have used it to successfully separate mitogenomes (~16,000 bp) of two shark species, whose reads (from gene capture, paired-end, read length 300bp) have been combined intentionally. Long-read sequencing is the better way to go for separating long gene copies.

# How to generate the input data for copyseparator
1. PCR amplicons --- prepare libraries --- NGS short-read sequencing
2. Target gene capture --- NGS short-read sequencing
3. Some other ways ...

# Why short reads need to be mapped to a reference before runing copyseparator?
Mapping to a reference can organize the short reads in order. In this way, the subsets, the picked clusters and their respective consensus sequences, and the assembled gene copy sequences can all be in the correct order.

# Why the read-length has to be 250bp or longer?
During assembling, short reads can easily result in chimeric sequences. Contact me if you really want to use copyseparator for NGS data with a read length < 250bp.

# How to assemble gene copies if "copy_assemble" fails
I use SeaView (http://doua.prabi.fr/software/seaview) to examine the results from "copy_separate". I prefer to assemble gene copies by moving sequences around and link gene copies by eye. It just like assembling the forward and reverse sequences from Sanger sequencing. Even if you prefer to use "copy_assemble", you still need to check the alignment first. Pay special attention to the nucleotide overhangs introduced by mistake during the calculation of consensus sequences of picked clusters.

# What if one or a few subsets did not get the anticipated number of gene copies?
If "copy_number" is set as 2, copyseparator will pick two of the largest clusters and assume that they represent the two gene copies. Most times, this is the case. Occasionally, for one or a few subsets, the picked two clusters both represent the same gene copy. When this happens, below is what you may do.

1. If the value of "overlap" is high (lots of overlaps between subsets), you may ignore these subsets and can still assemble the gene copies successfully.

2. You may set "copy_number" as 3 or 4 and re-run "copy_separate".

3. Or you can run "copy_detect" for specific subsets by setting the "copy_number" as 3 or 4.

# How to cite copyseparator
Run the following to get the appropriate citation for the version you’re using:

  citation(package = "copyseparator")

To cite package ‘copyseparator’ in publications use:

  Yang L (2022). _copyseparator: Assembling Long Gene Copies from Short Read Data_. R package version 1.0.0,
  <https://CRAN.R-project.org/package=copyseparator>.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {copyseparator: Assembling Long Gene Copies from Short Read Data},
    author = {Lei Yang},
    year = {2022},
    note = {R package version 1.0.0},
    url = {https://CRAN.R-project.org/package=copyseparator},
  }

ATTENTION: This citation information has been auto-generated from the package DESCRIPTION file and may need manual editing, see
‘help("citation")’.
