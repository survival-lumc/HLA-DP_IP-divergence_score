#' @title Create a patient-donor venn diagram
#' 
#' @param Venn_Data data frame with for each row a specific peptide and for each column an allele and whether or not this peptide is present in the allele
#' @param Pat1 string. An allele of the patient
#' @param Pat2 string. Second allele of the patient
#' @param Don1 string. An allele of the donor
#' @param Don2 string. Second allele of the donor
#' 
#' @import ggVennDiagram
#' 
#' @author Lars van der Burg
venn_diagram = function(sequence_data, Pat1, Pat2, Don1, Don2){
  
  if(FALSE %in% (c("Sequence", Pat1, Pat2, Don1, Don2) %in% colnames(sequence_data))){
    stop(glue("Please make sure that the columns Sequence, {Pat1}, {Pat2}, {Don1} and {Don2} are present in the sequence_data"))
  }
  
  
  Venn_data = sequence_data |> 
    mutate(Pat1 = ifelse(!!sym(Pat1) != 0, Sequence, NA),
           Pat2 = ifelse(!!sym(Pat2) != 0, Sequence, NA),
           Don1 = ifelse(!!sym(Don1) != 0, Sequence, NA),
           Don2 = ifelse(!!sym(Don2) != 0, Sequence, NA)) |> 
  select(Sequence, Pat1, Pat2, Don1, Don2)
    
  
  don_pat = list(Venn_data[["Don1"]], Venn_data[["Don2"]], Venn_data[["Pat1"]], Venn_data[["Pat2"]]); names(don_pat) = c(Don1, Don2, Pat1, Pat2)
  don_pat = lapply(don_pat, function(x){x[!is.na(x)]})
  
  A = Venn(don_pat); B = process_data(A, shape_id = "401f")

  
  colors <- c("1" = "#FF3333", "2" = "#FF3333", "3" = "#FFFF00", "4" = "#FFFF00", "1/2" = "#FF3333", "1/3" = "#0099CC", "1/4" = "#0099CC", "2/3" = "#0099CC",
              "2/4" = "#0099CC", "3/4" = "#FFFF00", "1/2/3" = "#0099CC", "1/2/4" = "#0099CC", "1/3/4" = "#0099CC", "2/3/4" = "#0099CC", "1/2/3/4" = "#0099CC")
  
  legend_colors <- c("Donor-specific peptides" = "#FFFF00", "Patient-specific peptides" = "#FF3333", "Shared peptides" = "#0099CC")
  legend_data <- data.frame(legend_id = names(legend_colors), x = c(0.5, 0.5, 0.5), y = c(0.5, 0.5, 0.5))
  
  
  p = ggplot() + 
    # This is just to get a legend, other layers override the points
    geom_point(data = legend_data, aes(x, y, fill = legend_id), shape = 22, size = 5, show.legend = TRUE) +
    
    geom_polygon(aes(X, Y, fill = id, group = id), data = venn_regionedge(B), show.legend = FALSE) + 
    
    geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(B), show.legend = FALSE, linewidth = 1, color = "black") +  
    geom_text(aes(X, Y, label = paste0(name, " (", str_replace_all(count, pattern = "(\\d{3}$)", replacement = paste0(",", unlist(str_extract_all(count, "(\\d{3}$)")))), ")")), data = venn_setlabel(B)) +
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

