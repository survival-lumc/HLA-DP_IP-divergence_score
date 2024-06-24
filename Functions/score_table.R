#' @title Give for a patient-donor combination all calculated scores
#' 
#' @param mismatches data.frame. Containing for all theoretical possible allele combinations the calculated scores.
#' @param Pat1 string. An allele of the patient
#' @param Pat2 string. Second allele of the patient
#' @param Don1 string. An allele of the donor
#' @param Don2 string. Second allele of the donor
#'
#' @return A table with for that combination the HLA-DPB1 mismatch scores.
#' 
#' 
#' @author Lars van der Burg
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

  tab = mismatches[mismatches$Donor == paste(Don1, Don2, sep = "+") & mismatches$Patient == paste(Pat1, Pat2, sep = "+"), ]
  
  return(tab)
}
