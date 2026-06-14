#' @title Give for a patient-donor combination all calculated scores
#' 
#' @param IP_scores data.frame. Containing for all theoretical possible allele combinations the calculated scores.
#' @param Pat1 string. An allele of the patient
#' @param Pat2 string. Second allele of the patient
#' @param Don1 string. An allele of the donor
#' @param Don2 string. Second allele of the donor
#'
#' @return A table with for that combination the HLA-DPB1 IP-divergence scores.
#' 
#' 
#' @author Lars van der Burg
score_table = function(IP_scores, Pat1, Pat2, Don1, Don2){
  
  ## Determine DT Patient
  DT_Pat = paste(Pat1, Pat2, sep = "+")
  index_Pat = IP_scores$Patient %in% DT_Pat
  
  ### Order is important, maybe other configuration
  if(sum(index_Pat) == 0){
    DT_Pat = paste(Pat2, Pat1, sep = "+")
    index_Pat = IP_scores$Patient %in% DT_Pat
  }
  
  
  ## Determine DT Donor
  DT_Don = paste(Don1, Don2, sep = "+")
  index_Don = IP_scores$Donor %in% DT_Don
  
  ### Order is important, maybe other configuration
  if(sum(index_Don) == 0){
    DT_Don = paste(Don2, Don1, sep = "+")
    index_Don = IP_scores$Donor %in% DT_Don
  }
  

  ## IP-score  
  IP_score = IP_scores[index_Pat & index_Don, ]

  
  if(nrow(IP_score) == 0){
    message("With the supplied combination of alleles there is no IP-divergence score available.")
  }

  
  return(IP_score)
}
