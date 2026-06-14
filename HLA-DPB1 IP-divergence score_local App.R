## Start -----
## Code to create the IP-divergence score Shiny app locally
## Download the required packages and run all code
#
## written by: Lars van der Burg
#

## Load required packages -----
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  
  library(tidyverse)
  library(ggpubr)
  library(ggVennDiagram)
  library(RColorBrewer)
  
  library(plotly)
})


## Function to filter the scores of a specific Patient-Donor combination
source("Functions/score_table.R")

## Function to create a Venn diagram based on the sequence data for a specific Patient-Donor combination
source("Functions/Venn_diagram.R")


## Read-in the sequence data
load("Data/sequence_data.RData")

## Read-in pre-calculated IP-divergence scores to circumvent long loading times in the Shiny app
load("Data/IP_scores.RData")



## Define UI for dataset viewer app -----
ui <- fluidPage(theme = shinytheme("readable"),
  
  titlePanel("HLA-DP IP-divergence scores"),
  
  helpText('The immunopeptidomes (IPs) were measured for K562 cell lines transduced with 15 distinct HLA-DP allotypes, with the data publicly available in the PRIDE partner repository under dataset identifier PXD030591. 
           HLA-DP IP-divergence scores were derived based on the difference between the IPs of patient and donor-specific allotype. 
           The "count score" ranges from 0-8520 and represents the total number of peptides specific to the patient, while the "ratio score" ranges from 0-1 and is calculated as the proportion of patient-specific peptides relative to the entire patient-derived peptide repertoire, thereby accounting for the heterozygosity of the patient. 
           Thereby, the HLA-DP IP-divergence scores quantify the magnitude of the immunogenicity of a HLA-DPB1 mismatch between patient and donor in the graft-versus-host direction.'),
  
  sidebarLayout(
    
    sidebarPanel(
      selectInput(inputId = "Don1",
                  label = "Donor HLA-DPB1 allele 1:",
                  choices = c("01:01", "02:01", "03:01", "04:01", "04:02", "05:01", "06:01", "09:01", "10:01", "11:01", "13:01", "14:01", "15:01", "17:01", "19:01"), 
                  selected = "01:01"),
      
      selectInput(inputId = "Don2",
                  label = "Donor HLA-DPB1 allele 2:",
                  choices = c("01:01", "02:01", "03:01", "04:01", "04:02", "05:01", "06:01", "09:01", "10:01", "11:01", "13:01", "14:01", "15:01", "17:01", "19:01"),
                  selected = "02:01"),

      selectInput(inputId = "Pat1",
                  label = "Patient HLA-DPB1 allele 1:",
                  choices = c("01:01", "02:01", "03:01", "04:01", "04:02", "05:01", "06:01", "09:01", "10:01", "11:01", "13:01", "14:01", "15:01", "17:01", "19:01"),
                  selected = "03:01"),

      selectInput(inputId = "Pat2",
                  label = "Patient HLA-DPB1 allele 2:",
                  choices = c("01:01", "02:01", "03:01", "04:01", "04:02", "05:01", "06:01", "09:01", "10:01", "11:01", "13:01", "14:01", "15:01", "17:01", "19:01"),
                  selected = "04:01")),
    
    mainPanel(
      textOutput("info"),
      tableOutput("table"),
      plotOutput("plot")
    )
  )
)


## Define server logic to summarize and view selected dataset -----
server <- function(input, output) {
  
  output$table = renderTable(score_table(IP_scores, input$Pat1, input$Pat2, input$Don1, input$Don2))
  output$plot = renderPlot(venn_diagram(sequence_data, input$Pat1, input$Pat2, input$Don1, input$Don2))
}


## Create Shiny app -----
shinyApp(ui = ui, server = server)