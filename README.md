# HLA-DPB1 immunopeptidome divergence score
R code repository for the HLA-DP immunopeptidome (IP-)divergence scores described in the manuscript ...


# Usage
There is an interactive shiny application of the HLA-DPB1 IP-divergence score tool online which can be used: https://jgkvanderhem.shinyapps.io/HLA-DP_IP-divergence_scores/

It is also possible to run the tool locally on your own laptop. If you are a git user, you can clone the repository by

`git clone https://github.com/survival-lumc/HLA-DP_IP-divergence_score.git`

Otherwise, you can simply download a zip file containing the directory by clicking Code -> Download ZIP at the top-right of this Github page. Extract the zipped files to a directory of your choice.

Afterwards, you can double-click the HLA-DP_IP-divergence_score.Rproj file to open an Rstudio session in the directory you have just downloaded. This will ensure all file-paths called in the files are maintained. The minimal script and the .Rmd files can now be executed.


# Contains
The repository contains the following:
* HLA-DPB1 IP-divergence score_tutorial.Rmd: script that shows how the data looks like, how the IP-divergence scores are calculated and what a Patient-Specific  to run the interactive shiny application locally. 
  * sequence_data.RData: sequence_data as described in the manuscript ...
  * IP_divergence_score.R: function that creates the IP-divergence scores based on the sequence data
  * Venn_diagram.R: function that creates the Venn diagram for a chosen Patient-Donor combination
* HLA-DPB1 IP-divergence score_local App.R: script to create the HLA-DPB1 IP-divergence tool locally.
