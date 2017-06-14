init_ext_valid_obj <- function(model_CV_obj, random = F, train_data_set = c("K-array", "H-array"), train_position = "X", test_position = "Z") {
  # this function initializes the ext_valid_obj
  # this object is used for external cross validation for the function master_ext_valid
  
  source("support_functions/load_AA_property_data.R")
  
  data_sets = c("K-array", "H-array")
  
  # load training data
  temp_list = load_AA_property_data(data_set = train_data_set, position = train_position, random = random)
  train_data = temp_list[[1]]
  train_response = temp_list[[2]]
  
  # load testing data
  temp_list = load_AA_property_data(data_set = data_sets[!(data_sets %in% train_data_set)], position = test_position, random = random)
  test_data = temp_list[[1]]
  test_response = temp_list[[2]]
  
  # remove the colnames because it won't make sense with the cross validation
  colnames(train_data) = 1:ncol(train_data)
  colnames(test_data) = 1:ncol(test_data)
  
  ###########################
  # stitch into a list
  
  ext_valid_obj = list()
  ext_valid_obj$train_data = train_data
  ext_valid_obj$test_data = test_data
  ext_valid_obj$train_response = train_response
  ext_valid_obj$test_response = test_response
  ext_valid_obj$model_func = model_CV_obj$model_func
  ext_valid_obj$predict_func = model_CV_obj$predict_func
  ext_valid_obj$random_response_func = model_CV_obj$random_response_func
  
  return(ext_valid_obj)
  
}
