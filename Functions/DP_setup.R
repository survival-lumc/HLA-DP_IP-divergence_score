#' @title Function to calculate mismatch scores
#'
#' @description Calculate the two mismatch scores (count and ratio) for the supplied genotyings. A mismatch score can only be calculated if the peptides of all four alleles are known
#'
#' @param Data data.frame containing the HLA-DP information
#' @param DPB1_genes vector indicating the four genenames (two patient, two donor) as they are called in `Data`
#' @param peptide_data the available peptide data for all HLA alleles
#' 
#' @return A list with two elements: 1) the calculated mismatch scores and 2) an index for which individuals no mismatch score could be calculated
#' 
#' 
#' @author Lars van der Burg
DP_setup = function(Data, DPB1_genes = c("HLA_DPP1", "HLA_DPP2", "HLA_DPD1", "HLA_DPD2"), peptide_data){
  
  N = nrow(Data)

## Retrieve HLA-DP data from the data    
  haplos = unique(sort(as.matrix(Data[, DPB1_genes])))
  haplos = haplos[order(unlist(lapply(strsplit(haplos, ":"), function(x){as.numeric(str_sub(x[2], end = 2))})))]; 
  haplos = haplos[order(unlist(lapply(strsplit(haplos, ":"), function(x){as.numeric(x[1])})))]

  
  DTs_don = apply(Data[, DPB1_genes[3:4]], 1, function(x){paste(sort(x), collapse = "+")})
  DTs_pat = apply(Data[, DPB1_genes[1:2]], 1, function(x){paste(sort(x), collapse = "+")})
  DTs_match = apply(cbind(DTs_don, DTs_pat), 1, function(x){paste(x, collapse = "&")})
  
  
## Put peptide data in long format
  peptide_data = cbind("Allele1" = colnames(peptide_data), peptide_data)  # Allele1 = donor
  peptide_data_long = gather(peptide_data, "Allele2", "Mismatch", -Allele1)  # Allele2 = patient
  
  peptide_data_long = peptide_data_long[order(peptide_data_long$Allele2), ]
  peptide_data_long = peptide_data_long[order(peptide_data_long$Allele1), ]
  rownames(peptide_data_long) = 1:nrow(peptide_data_long)
  
  
## Calculate mismatch score for data
  MMS = rep(NA, N); discard_index = NULL
  for(i in 1:N){
    cat("Now determining the score for individual:", i, "of the", N, "\r")
  
    alleles_don = DTs_don[[i]]; splt_alleles_don = unlist(strsplit(alleles_don, "\\+")) 
    alleles_pat = DTs_pat[[i]]; splt_alleles_pat = unlist(strsplit(alleles_pat, "\\+"))
      
    index_allele1 = which(peptide_data_long$Allele1 == alleles_don)
    index_allele2 = which(peptide_data_long$Allele2 == alleles_pat)
      
    if(length(index_allele1) != 0 & length(index_allele2) != 0){
      MMS[i] = peptide_data_long$Mismatch[index_allele1[index_allele1 %in% index_allele2]]
    } else {
      discard_index = c(discard_index, i)
    }
  }  
    

  return(list(mismatch_score = MMS, discard_index = discard_index))
}
