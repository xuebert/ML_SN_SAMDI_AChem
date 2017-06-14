# production script for sampling peptides and predicting SN of unsampled peptides (Figure 6)
num_boot = 1 # number of bootstraps.  Need to parallelize script (or wait a long time) to run many bootstraps

source("ML_functions/init_model_CV_obj.R")
source("ML_functions/init_ext_valid_obj.R")
source("ML_functions/master_ext_valid.R")

# specify datasets
data_sets = list("K-array", "H-array")

# initialize storage variables
min_samples = 2
storage_template = matrix(NA, nrow = num_boot, ncol = 361) # columns correspond to number of sampled peptides.  Rows are for each bootstrap
storage_template = list(K = storage_template, H = storage_template) # for each K- and H-array
storage_template = list(actual = storage_template, random = storage_template)
Q2_sampled = storage_template

# storage information
# top list layer holds actual and random (control) data
# second layer holds K and H arrays
# lowest layer has matrices where rows are bootstraps and columns are number of sampled peptides

for (n_data_set in 1:length(data_sets)) {
  print("#######################")
  print(paste("dataset", data_sets[[n_data_set]]))
  
  # initialize objects
  model_CV_obj = init_model_CV_obj(random = F, data_set = data_sets[[n_data_set]])
  ext_valid_obj = init_ext_valid_obj(model_CV_obj, random = F, train_data_set = data_sets[[n_data_set]], train_position = c("X", "Z"), test_position = c("X", "Z"))
  
  # format train data.  "external" dataset is the non-sampled data
  # store the full datasets
  train_data = ext_valid_obj$train_data
  train_response = ext_valid_obj$train_response
  
  # iterate over sample sizes
  for (n_samples in 360:min_samples) {
    print(paste("sample ", n_samples, sep = ""))
    
    # re-randomize for each sample
    model_CV_obj_random = init_model_CV_obj(random = F, data_set = data_sets[[n_data_set]])
    ext_valid_obj_random = init_ext_valid_obj(model_CV_obj_random, random = T, train_data_set = data_sets[[n_data_set]], train_position = c("X", "Z"), test_position = c("X", "Z"))
    
    # format train data.  "external" dataset is the non-sampled data
    # store the full datasets
    train_data_random = ext_valid_obj_random$train_data
    train_response_random = ext_valid_obj_random$train_response
    
    # iterate over the specified number of bootstraps
    for (n_boot in 1:num_boot) {
      
      # select train/test samples
      select_samples = sample(1:361, n_samples)
      ext_valid_obj$train_data = train_data[select_samples, ]
      ext_valid_obj$train_response = train_response[select_samples]
      ext_valid_obj$test_data = train_data[-select_samples, ]
      ext_valid_obj$test_response = train_response[-select_samples]
      
      ext_valid_obj_random$train_data = train_data_random[select_samples, ]
      ext_valid_obj_random$train_response = train_response_random[select_samples]
      ext_valid_obj_random$test_data = train_data_random[-select_samples, ]
      ext_valid_obj_random$test_response = train_response_random[-select_samples]
      
      # calculate Q2 
      Q2_sampled$actual[[n_data_set]][[n_boot, n_samples]] = master_ext_valid(ext_valid_obj)[[1]]
      Q2_sampled$random[[n_data_set]][[n_boot, n_samples]] = master_ext_valid(ext_valid_obj_random)[[1]]
      
    }
  }
}
filename = paste("sampling_peptides_results.RData", sep = "")

# NOTE: original sampling_peptides_results.RData was run on Northwestern's HPC, Quest.  This RData file is in a different format as a result.  Instead of storing all bootstraps, confidence intervals were calculated and stored.  If you run your own code, your result will store each individual bootstrap.  Your RData file then will not work properly with make_figure6.R
save(Q2_sampled, file = filename)
