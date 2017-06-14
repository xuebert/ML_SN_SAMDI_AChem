# this is the top level production_script for predicting signal to noise ratio using random forests
# this does three datasets: K-array (KXZ), H-array, and the combination of both
# this generates the results behind Figure 5 and STable 1 and 2

rm(list=ls())

source("ML_functions/init_model_CV_obj.R") # sets up internal model cross validation object
source("ML_functions/master_CV.R") # runs internal k-folds CV
source("ML_functions/init_ext_valid_obj.R") # sets up external validation object
source("ML_functions/master_ext_valid.R") # runs external validation

num_AAPP = 39
# select feature sets
PP_lists = list(1:num_AAPP, (num_AAPP+1):(2*num_AAPP), 1:(2*num_AAPP)) # x-position all AAPP, x-position all AAPP, both positions all AAPP
PP_lists = c(PP_lists, 1:(2*num_AAPP), lapply(1:num_AAPP, function(i) c(i, i+num_AAPP)))# each individual AAPP in x position, z position, and both positions

# specify datasets
data_sets = list("K-array", "H-array", c("K-array", "H-array"))

# initialize storage variables
storage_template = vector("list", 6)
dim(storage_template) = c(3,2) # 1st/2nd row are individual datasets.  3rd row is both datasets.  1st/2nd columns are all physical properties, then individual physical properties
Q2_int = storage_template
Q2_int_random = storage_template
Q2_ext = matrix(NA, ncol = 2, nrow = 5, dimnames = list(c("trainX_testX", "trainX_testZ", "trainZ_testX", "trainZ_testZ", "trainXZ_testXZ"), c("K-array_train", "H-array_train")))
Q2_ext_random = Q2_ext

# Q2_int holds internal k-folds results within arrays.  List rows are K-array, H-array, and both arrays.
# list columns are (1) all features at once and (2) individual features
# first list column has three values: x-position, z-position, and both positions
# second list column has 1:39 for x-position, 40:78 for x-position, and 79:117 for both positions

# specify train/test position for external Q2
train_positions = list("X", "X", "Z", "Z", c("X", "Z"))
test_positions = list("X", "Z", "X", "Z", c("X", "Z"))

for (n_data_set in 1:length(data_sets)) {
  
  print("################################")
  print(paste("dataset ", data_sets[[n_data_set]]))
  
  #### internal Q2 for separate datasets.  Needed to perform k-folds internal validation
  # also needed to setup external validation
  # initialize objects
  model_CV_obj = init_model_CV_obj(random = F, data_set = data_sets[[n_data_set]])
  model_CV_obj_random = init_model_CV_obj(random = T, data_set = data_sets[[n_data_set]])
  
  # keep track of data and response
  data_mat = model_CV_obj$data_mat
  data_mat_random = model_CV_obj_random$data_mat
  
  # calculate Q2 for all AA physical properties.  This function takes a few seconds to run
  for (n_PP in 1:3) {
    model_CV_obj$data_mat = data_mat[,PP_lists[[n_PP]]]
    Q2_int[[n_data_set, 1]] = c(Q2_int[[n_data_set, 1]], mean(master_CV(model_CV_obj)[[1]]))

    model_CV_obj_random$data_mat = data_mat_random[,PP_lists[[n_PP]]]
    Q2_int_random[[n_data_set, 1]] = c(Q2_int_random[[n_data_set, 1]], mean(master_CV(model_CV_obj_random)[[1]]))
  }
  
  # calculate Q2 for individual AA physical properties.  This function takes quite a few minutes to run.  Was originally integrated and run on Northwestern's high performance computing cluster (Quest)
  for (n_PP in 4:length(PP_lists)) {
    model_CV_obj$data_mat = data_mat[,PP_lists[[n_PP]], drop = F]
    Q2_int[[n_data_set, 2]] = c(Q2_int[[n_data_set, 2]], mean(master_CV(model_CV_obj)[[1]]))

    model_CV_obj_random$data_mat = data_mat_random[,PP_lists[[n_PP]], drop = F]
    Q2_int_random[[n_data_set, 2]] = c(Q2_int_random[[n_data_set, 2]], mean(master_CV(model_CV_obj_random)[[1]]))
  }
  
  #### external Q2 (training on one array, testing on the other)
  if (n_data_set != 3) { # ignore combined dataset for external validation
    
    # iterate over different combinations of positions
    for (n_external in 1:5) {
      # do the external CV
      ext_valid_obj = init_ext_valid_obj(model_CV_obj, random = F, train_data_set = data_sets[[n_data_set]], train_position = train_positions[[n_data_set]], test_position = test_positions[[n_data_set]])
      ext_valid_obj_random = init_ext_valid_obj(model_CV_obj_random, random = T, train_data_set = data_sets[[n_data_set]], train_position = train_positions[[n_data_set]], test_position = test_positions[[n_data_set]])
      
      # run the external validation and store external Q2
      Q2_ext[[n_external, n_data_set]] = master_ext_valid(ext_valid_obj)
      Q2_ext_random[[n_external, n_data_set]] = master_ext_valid(ext_valid_obj_random)
    }
  }
}

save(Q2_ext, Q2_int, Q2_ext_random, Q2_int_random, file = "predict_SN_results.RData")
