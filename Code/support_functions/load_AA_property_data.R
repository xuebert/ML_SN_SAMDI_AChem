load_AA_property_data <- function(data_set = c("K-array", "H-array"), position = "X", random = F) {
  # this function loads the AA_physical_properties data for K-array or H-array in either X, Z or both positions
  
  source("../Data/amino_acid_physical_properties/AA_physical_properties.R")
  
  # storage variables
  response = c()
  peptides = c()
  data_mat = matrix(NA, ncol = 0, nrow = 361)
  
  if ("K-array" %in% data_set) {
    
    load("../Data/K_array_data.RData")
    SN_K = colMeans(SN_K)
    response = c(response, SN_K)
    
    AA1 = rep(c("A", "R", "N", "D", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 19)
    AA2 = rep(c("A", "R", "N", "D", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), each = 19)
    peptides = paste(AA1, AA2, sep="")
    
    # add in the AA_physical_properties
    data_mat = AA_physical_properties(AA1)
    colnames(data_mat) = paste("K_", "X_", colnames(data_mat), sep="")
    data_mat_temp = AA_physical_properties(AA2)
    colnames(data_mat_temp) = paste("K_", "Z_", colnames(data_mat_temp), sep="")
    data_mat = cbind(data_mat, data_mat_temp)
    rownames(data_mat) = peptides
    
    # retain the proper positions
    if (!("X" %in% position)) { # if X is not in train positions
      # keep second half of columns
      data_mat = data_mat[, (ncol(data_mat)/2 + 1):ncol(data_mat)]
    } else if (!("Z" %in% position)) { # if Z is not in train positions
      # keep first half of columns
      data_mat = data_mat[, 1:(ncol(data_mat)/2)]
    }

  }
  
  if ("H-array" %in% data_set) {
    
    load("../Data/H_array_data.RData")
    SN_H = colMeans(SN_H)
    response = c(response, SN_H)
    
    AA1 = rep(c("A", "R", "N", "D", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), 19)
    AA2 = rep(c("A", "R", "N", "D", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), each = 19)
    peptides = paste(AA1, AA2, sep="")
    
    # add in the AA_physical_properties
    data_mat2 = AA_physical_properties(AA1)
    colnames(data_mat2) = paste("H_", "X_", colnames(data_mat2), sep="")
    data_mat2_temp = AA_physical_properties(AA2)
    colnames(data_mat2_temp) = paste("H_", "Z_", colnames(data_mat2_temp), sep="")
    data_mat2 = cbind(data_mat2, data_mat2_temp)
    rownames(data_mat2) = peptides
    
    # retain the proper positions
    if (!("X" %in% position)) { # if X is not in train positions
      # keep second half of columns
      data_mat2 = data_mat2[, (ncol(data_mat2)/2 + 1):ncol(data_mat2)]
    } else if (!("Z" %in% position)) { # if Z is not in train positions
      # keep first half of columns
      data_mat2 = data_mat2[, 1:(ncol(data_mat2)/2)]
    }
    
    # combine with phosphatase if there
    if ("K-array" %in% data_set) {
      data_mat = rbind(data_mat, data_mat2)
    } else {
      data_mat = data_mat2
    }
    
  }
  
  # randomize if specified
  if (random) {
    
    num_rows = nrow(data_mat)
    
    for (n_col in 1:ncol(data_mat)) {
      data_mat[, n_col] = data_mat[sample(1:num_rows, num_rows), n_col]
    }
  }
  
  return(list(data_mat, response))
}
