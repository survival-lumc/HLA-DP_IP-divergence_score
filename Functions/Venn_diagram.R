#' @title Create a patient-donor venn diagram
#' 
#' @param Venn_Data data frame with for each row a specific peptide and for each column an allele and whether or not this peptide is present in the allele
#' @param Pat1 string. An allele of the patient
#' @param Pat2 string. Second allele of the patient
#' @param Don1 string. An allele of the donor
#' @param Don2 string. Second allele of the donor
#' 
#' 
#' @author Lars van der Burg
venn_diagram = function(Venn_Data, Pat1, Pat2, Don1, Don2){
  
  pat1 = paste("DP", paste(unlist(strsplit(as.character(Pat1), ":")), collapse = ""), sep = ""); pat2 = paste("DP", paste(unlist(strsplit(Pat2, ":")), collapse = ""), sep = "")
  don1 = paste("DP", paste(unlist(strsplit(Don1, ":")), collapse = ""), sep = ""); don2 = paste("DP", paste(unlist(strsplit(Don2, ":")), collapse = ""), sep = "")
  
  don_pat = list(Venn_Data[, don1], Venn_Data[, don2], Venn_Data[, pat1], Venn_Data[, pat2]); names(don_pat) = c(Don1, Don2, Pat1, Pat2)
  A = Venn(don_pat); B = process_data(A, shape_id = "401f")

  venn_numbers = venn_regionlabel(B); venn_numbers[15, "count"] = venn_numbers[15, "count"] - 1
  
  
  colors <- c("1" = "#7570B3", "2" = "#7570B3", "3" = "#D95F02", "4" = "#D95F02", "1/2" = "#7570B3", "1/3" = "#1B9E77", "1/4" = "#1B9E77", "2/3" = "#1B9E77",
              "2/4" = "#1B9E77", "3/4" = "#D95F02", "1/2/3" = "#1B9E77", "1/2/4" = "#1B9E77", "1/3/4" = "#1B9E77", "2/3/4" = "#1B9E77", "1/2/3/4" = "#1B9E77")
  
  
  p = ggplot() + 
    geom_polygon(aes(X, Y, fill = id, group = id), data = venn_regionedge(B), show.legend = FALSE) + 
    
    geom_path(aes(X, Y, color = id, group = id), data = venn_setedge(B), show.legend = FALSE, linewidth = 1, color = "black") +  
    geom_text(aes(X, Y, label = name), data = venn_setlabel(B)) +
    geom_label(aes(X, Y, label = count), data = venn_numbers) +
    
    scale_fill_manual(values = colors) + 
    coord_equal() +
    
    annotate("text", label = "Donor", x = 0.125, y = 0.85, size = 5, fontface = 2) +
    annotate("text", label = "Patient", x = 0.875, y = 0.85, size = 5, fontface = 2) +
    
    theme_void()
  
  return(p)
}

