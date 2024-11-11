# HLA-DPB1 immunopeptidome mismatch score
R code repository for the HLA-DPB1 immunopeptidome mismatch score described in the manuscript ...

The repository contains the following:
* Mismatch Score tutorial.Rmd: a simple tutorial showing how to obtain the mismatch scores and Venn diagram for a specific patient-donor combinations. By adjusting the specific alleles the patient or donor has, the different scores can be obtained
* Peptide_data_count.Rdata and Peptide_data_ratio.Rdata: datafile containing the HLA-DPB1 immunopeptidome count and ratio mismatch scores for all measured patient-donor combinations
* Venn_data.Rdata: the different peptides that are measured for every HLA-DPB1 allele. Every row is a unique peptide which is measured for the allele if it is not NA.
* Additional code and an example dataset used to calculate create the scores, make the Venn diagram or allow for more insights in the data.  

# Usage
There is an interactive shiny application of the HLA-DPB1 mismatch tool online which can be used: https://lljvdburg.shinyapps.io/myapp/

It is also possible to run the tool locally on your own laptop. If you are a git user, you can clone the repository by

`git clone https://github.com/survival-lumc/HLA-DP_Mismatch_scores.git`

Otherwise, you can simply download a zip file containing the directory by clicking Code -> Download ZIP at the top-right of this Github page. Extract the zipped files to a directory of your choice.

Afterwards, you can double-click the HLA-DP_Mismatch_scores.Rproj file to open an Rstudio session in the directory you have just downloaded. This will ensure all file-paths called in the files are maintained. The minimal script and the .Rmd files can now be executed.


