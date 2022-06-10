# How to install copyseparator
install.packages("copyseparator")
  
# Examples on how to run copyseparator
https://github.com/LeiYang-Fish/copyseparator_more_examples

# What if one or a few subsets did not get the anticipated number of gene copies?
If "copy_number" is set as 2, copyseparator will pick two of the largest clusters and assume that they represent the two gene copies. Most times, this is the case. Occasionally, for one or a few subsets, the picked two clusters both represent the same gene copy.

You may set "copy_number" as 3 or 4 and re-run "copy_separate".

Or you can run "copy_detect" for specific subsets by setting the "copy_number" as 3 or 4.

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
