AA_physical_properties <- function(AA_list, full_name = T) {
  # this function returns the physical properties for the specified amino acids.  AA_list can be a list or char vector
  # from this source: http://onlinelibrary.wiley.com/store/10.1002/bip.20296/asset/supinfo/jwsBIP.mei.pdf?v=1&s=3acfbe1500c355aeb4a9a36a25fafe0d5e676c37
  # amino acids need to be specified as capital letters
  
  # load amino acid physical properties
  load("../Data/amino_acid_physical_properties/AA_physical_properties.RData")
  AA_names = rownames(AA_physical_properties_data)
  
  # storage
  AA_physical_properties_data_return = matrix(NA, ncol = 50, nrow = length(AA_list))
  # iterate and load physical properties for each specified amino acid
  for (n_AA in 1:length(AA_list)) {
    AA_physical_properties_data_return[n_AA, ] = AA_physical_properties_data[which(AA_names %in% AA_list[[n_AA]]), ]
  }
  
  # load property names
  property_names = read.table("../Data/amino_acid_physical_properties/AA_property_names.txt", header = F, sep = "!", colClasses = "character")[[1]]
  # store property names
  if (!full_name) {
    property_names = 1:length(property_names)
  }
  colnames(AA_physical_properties_data_return) <- property_names
  
  # remove uninformative ones
  removal = c(1, 6, 7, 12, 13:18, 34)
  AA_physical_properties_data_return = AA_physical_properties_data_return[, -removal]
  
  return(AA_physical_properties_data_return)
  
}
