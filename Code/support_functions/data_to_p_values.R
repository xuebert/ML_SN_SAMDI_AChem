data_to_p_values <- function(data_value, AA_storage, low_threshold = 30, high_threshold = 30, test = "Fischer") {
  # this function calculates p-values for AA in low and high regions
  # this function assumes that AA_storage is a matrix of amino acid letters and that data_value is the corresponding value (signal to noise) of that amino acid pair
  # the number of columns in AA_storage corresponds to the number of amino acid positions
  
  # make sure it is a matrix
  AA_storage = as.matrix(AA_storage)
  
  # sort by data_value
  data_order = order(data_value, decreasing = F)
  data_value = data_value[data_order]
  AA_storage = AA_storage[data_order, ]
  
  # the amino acids found across the arrays
  all_AAs = sort(unique(as.vector(AA_storage)))
  
  # record number of peptides and num of positions
  num_pep = length(data_value)
  num_positions = ncol(AA_storage)
  
  # fill out storage variables
  storage_table = matrix(0, nrow = 2 * num_positions, ncol = length(all_AAs), dimnames = list(paste(rep(c("low", "high"), each = num_positions), rep(1:num_positions, 2), sep=""), all_AAs)) # each row is a position to store AA information from the high or low regions
  AA_table = storage_table
  AA_table_p = storage_table
  AA_tot = matrix(0, ncol = length(all_AAs), nrow = num_positions)
  colnames(AA_tot) = all_AAs
  
  # store occurence of AAs within thresholds
  indices_low = 1:low_threshold # low threshold
  indices_high = num_pep:(num_pep - high_threshold) # high threshold
  
  ################# Fischer test ########################
  calculate_p <- function(num_pep, window_size, observed_in_window, peptides_of_interest) {
    return((sum(dhyper(observed_in_window:num_pep, peptides_of_interest, num_pep - peptides_of_interest, window_size))))
  }
  
  for (n_position in 1:num_positions) {
    
    # tabulate the AAs in the low region
    temp = table(AA_storage[indices_low, n_position])
    AA_table[n_position, names(temp)] = temp
    
    # same for high region
    temp = table(AA_storage[indices_high, n_position])
    AA_table[n_position + num_positions, names(temp)] = temp
    
    # count up total occurence of AA in n_position
    # first table all AAs
    temp = table(AA_storage[, n_position])
    # then introduce into storage
    AA_tot[n_position, names(temp)] = temp # this allows storage when AAs are missing
    
    ### calculate p values ###
    # low
    AA_table_p[n_position,] = sapply(1:length(all_AAs), function(n) calculate_p(num_pep=num_pep, window_size=low_threshold, observed_in_window=AA_table[[n_position,n]], peptides_of_interest=AA_tot[[n_position, n]]))
    
    # high
    AA_table_p[n_position + num_positions,] = sapply(1:length(all_AAs), function(n) calculate_p(num_pep=num_pep, window_size=high_threshold, observed_in_window=AA_table[[n_position + num_positions,n]], peptides_of_interest=AA_tot[[n_position, n]]))
    
  }
  
  return(list(AA_table_p, AA_table, AA_tot))
}