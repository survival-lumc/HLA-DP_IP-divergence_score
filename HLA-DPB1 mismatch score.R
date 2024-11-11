library(shiny)
library(shinythemes)

library(tidyverse)
library(ggpubr)
library(ggVennDiagram)
library(RColorBrewer)

library(plotly)


load("Data/Venn_Data.Rdata")
load("Data/Peptide_data_count.Rdata")
load("Data/Peptide_data_ratio.Rdata")


long_count = peptide_data_count |> 
  dplyr::mutate(Donor = rownames(peptide_data_count)) |> 
  tidyr::pivot_longer(-Donor) |> 
  dplyr::rename(Patient = name, Count = value)

long_ratio = peptide_data_ratio |> 
  dplyr::mutate(Donor = rownames(peptide_data_ratio)) |> 
  tidyr::pivot_longer(-Donor) |> 
  dplyr::rename(Patient = name, Ratio = value)

mismatches = long_count |> 
  dplyr::full_join(long_ratio, by = c("Donor", "Patient")) |> 
  dplyr::relocate(Patient) |> 
  dplyr::mutate(Count = as.integer(Count))



venn_diagram = function(Venn_Data, Pat1, Pat2, Don1, Don2){
  pat1 = paste("DP", paste(unlist(strsplit(as.character(Pat1), ":")), collapse = ""), sep = ""); pat2 = paste("DP", paste(unlist(strsplit(Pat2, ":")), collapse = ""), sep = "")
  don1 = paste("DP", paste(unlist(strsplit(Don1, ":")), collapse = ""), sep = ""); don2 = paste("DP", paste(unlist(strsplit(Don2, ":")), collapse = ""), sep = "")
  
  don_pat = list(Venn_Data[, don1], Venn_Data[, don2], Venn_Data[, pat1], Venn_Data[, pat2]); names(don_pat) = c(Don1, Don2, Pat1, Pat2)
  don_pat = lapply(don_pat, function(x){x[!is.na(x)]})
  
  A = Venn(don_pat); B = process_data(A, shape_id = "401f")

  
  colors <- c("1" = "#FFFF00", "2" = "#FFFF00", "3" = "#FF3333", "4" = "#FF3333", "1/2" = "#FFFF00", "1/3" = "#0099CC", "1/4" = "#0099CC", "2/3" = "#0099CC",
              "2/4" = "#0099CC", "3/4" = "#FF3333", "1/2/3" = "#0099CC", "1/2/4" = "#0099CC", "1/3/4" = "#0099CC", "2/3/4" = "#0099CC", "1/2/3/4" = "#0099CC")
  
  legend_colors <- c("Donor-specific peptides" = "#FFFF00", "Patient-specific peptides" = "#FF3333", "Shared peptides" = "#0099CC")
  legend_data <- data.frame(legend_id = names(legend_colors), x = c(0.5, 0.5, 0.5), y = c(0.5, 0.5, 0.5))
  
  
  p = ggplot() +
    # This is just to get a legend, other layers override the points
    geom_point(data = legend_data, aes(x, y, fill = legend_id), shape = 22, size = 5, show.legend = TRUE) +
    
    geom_polygon(aes(X, Y, fill = id, group = id), data = venn_regionedge(B), show.legend = FALSE) +

    geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(B), linewidth = 1, color = "black") +
    geom_text(aes(X, Y, label = paste0(name, "\n(", str_replace_all(count, pattern = "(\\d{3}$)", replacement = paste0(",", unlist(str_extract_all(count, "(\\d{3}$)")))), ")")), data = venn_setlabel(B)) +
    geom_label(aes(X, Y, label = count), data = venn_regionlabel(B)) +

    scale_fill_manual(values = c(colors, legend_colors), 
                      breaks = names(legend_colors), 
                      labels = names(legend_colors),
                      guide = guide_legend(title = "", override.aes = list(shape = 22, size = 5))) +

    coord_equal() +
    
    annotate("text", label = "Donor", x = 0.125, y = 0.85, size = 5, fontface = 2) +
    annotate("text", label = "Patient", x = 0.875, y = 0.85, size = 5, fontface = 2) +
    
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12))
    
  return(p)
}
score_table = function(mismatches, Pat1, Pat2, Don1, Don2){
  Pat01 = as.numeric(unlist(strsplit(Pat1, ":")))
  Pat02 = as.numeric(unlist(strsplit(Pat2, ":")))
  
  if((Pat01[1] > Pat02[1]) | (Pat01[1] == Pat02[1] && Pat01[2] > Pat02[2])){
    Pat_change = Pat1; Pat1 = Pat2; Pat2 = Pat_change
  } 
  
  
  Don01 = as.numeric(unlist(strsplit(Don1, ":")))
  Don02 = as.numeric(unlist(strsplit(Don2, ":")))
  
  if((Don01[1] > Don02[1]) | (Don01[1] == Don02[1] && Don01[2] > Don02[2])){
    Don_change = Don1; Don1 = Don2; Don2 = Don_change
  } 
  
  tab = mismatches[mismatches$Donor == paste(Don1, Don2, sep = "+") & mismatches$Patient == paste(Pat1, Pat2, sep = "+"), c("Count", "Ratio")]
  
  return(tab)
}



# Define UI for dataset viewer app ----
ui <- fluidPage(theme = shinytheme("readable"),
  
  titlePanel("HLA-DPB1 immunopeptidome mismatch score"),
  
  helpText('The immunepeptidomes (IPs) were measured for K562 cell lines transduced with 15 distinct HLA-DP allotypes, with the data publicly available in the PRIDE partner repository under dataset identifier PXD030591. 
           HLA-DPB1 immunopeptidome mismatch scores were derived based on the difference between the IPs of patient and donor-specific allotype. 
           The "count score" ranges from 0-8520) and represents the total number of peptides specific to the patient, while the "ratio score" ranges from 0-1 and is calculated as the proportion of patient-specific peptides relative to the entire patient-derived peptide repertoire, thereby accounting for the heterozygosity of the patient. 
           Thereby, the HLA-DPB1 immunopeptidome mismatch scores quantify the degree of HLA-DPB1 mismatch between patient and donor in the graft-versus-host direction.'),
  
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


# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  output$table = renderTable(score_table(mismatches, input$Pat1, input$Pat2, input$Don1, input$Don2))
  output$plot = renderPlot(venn_diagram(Venn_Data, input$Pat1, input$Pat2, input$Don1, input$Don2))
}


# Create Shiny app ----
shinyApp(ui = ui, server = server)