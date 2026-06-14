#' @title Function to calculate the IP-divergence scores 
#' 
#' @param sequence_data data.frame. The dataset with the measured peptides
#' 
#' @return A list with the following three elements:
#' - IP_score_count: Symmetrixc matrix of the count score
#' - IP_score_ratio: Symmetrixc matrix of the ratio score
#' - Overview: matrix with the patient and donor DTs, and the accompanying count and ratio scores
#' - crosstable_count: 
#' - crosstable_ratio: 
#' 
#' @author Lars van der Burg
#' 
IP_divergence_score = function(sequence_data){
  
# Settings
  ## Selection of the haplotypes
  HTs = sequence_data |>  
    dplyr::summarise_all(class) |> 
    tidyr::gather(allele, class) |> 
    
    filter(class %in% c("numeric", "integer")) |> 
    pull(allele)
  
  len_HTs = length(HTs)
  
  message("We will calculate the IP-divergence score for the following ", len_HTs, " alleles:\n", paste0(HTs, collapse = ", "), "\n")
  
  
  ## Convert possible NAs to 0, for calculation later
  sequence_data = sequence_data |> 
    mutate(across(all_of(HTs), \(x){ifelse(is.na(x), 0, x)}))
  
  
  ## Selection of the diplotypes
  DTs = apply(expand.grid(HTs, HTs), 1, function(x){paste(sort(x), collapse = "+")}) |> unique()
  len_DTs = length(DTs)
  splt_DTs = strsplit(DTs, "\\+")
  
  
  ## Selection of the peptides
  peptide_data = sequence_data |> select(all_of(HTs))
  
  nr_peptides = nrow(peptide_data)
  sum_peptides = colSums(peptide_data, na.rm = TRUE)
  


# Haplotype data
  cSums_HT = colSums(peptide_data, na.rm = TRUE)
  
  ## Scoring
  crosstable_count = crosstable_ratio = matrix(0, nrow = len_HTs, ncol = len_HTs, dimnames = list(HTs, HTs))
  for(i in 1:len_HTs){
    for(j in 1:len_HTs){
      
      if(i == j){
        crosstable_count[i, j] = sum(peptide_data[, i] != 0) 
      } else {
        crosstable_count[i, j] = sum(peptide_data[, i] != 0 & peptide_data[, j] != 0)
      }
      
      crosstable_ratio[i, j] = crosstable_count[i, j] / cSums_HT[i]
    }
  }
  
    
  
# Diplotype data
  peptide_data_DT = do.call(cbind, lapply(splt_DTs, \(x){apply(peptide_data[, x], 1, \(y){ifelse(TRUE %in% (y != 0), 1, 0)})}))
  colnames(peptide_data_DT) = DTs
  
  cSums_DT = colSums(peptide_data_DT, na.rm = TRUE)
  
  
  ## Scoring
  IP_score_count = IP_score_ratio = matrix(0, nrow = len_DTs, ncol = len_DTs, dimnames = list(DTs, DTs))
  for(i in 1:len_DTs){
    for(j in 1:len_DTs){
      
      ## Peptides in patient (j) not in the donor (i) are the problem
      if(i != j){
        IP_score_count[i, j] = sum(peptide_data_DT[peptide_data_DT[, i] == 0, j], na.rm = TRUE)  
        
        IP_score_ratio[i, j] = IP_score_count[i, j] / cSums_DT[j]
      }
    }
  }
  
  
  ## Reformating
  Count = IP_score_count |> as_tibble() |> 
    dplyr::mutate(Donor = rownames(IP_score_count)) |> 
    tidyr::pivot_longer(-Donor) |> 
    dplyr::rename(Patient = name, Count = value) |> 
    dplyr::mutate(Count = as.integer(Count))
  
  
  Ratio =  IP_score_ratio |> as_tibble() |>  
    dplyr::mutate(Donor = rownames(IP_score_ratio)) |> 
    tidyr::pivot_longer(-Donor) |> 
    dplyr::rename(Patient = name, Ratio = value)
  
  
  IP_scores = Count |> 
    dplyr::full_join(Ratio, by = c("Donor", "Patient")) |> 
    dplyr::relocate(Patient)
  
  
  
  return(list(IP_score_count = IP_score_count, IP_score_ratio = IP_score_ratio, 
              IP_scores = IP_scores, 
              crosstable_count = crosstable_count, crosstable_ratio = crosstable_ratio))
}
