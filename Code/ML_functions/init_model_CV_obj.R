init_model_CV_obj <- function(random = F, data_set = c("K-array", "H-array")) {
  # this function initializes the model_CV_obj
  # this object is used to do internal cross validation.  It is also required for setting up the external validation object.
  
  # ensure the right data_set names 
  if (!("K-array" %in% data_set) && !("H-array" %in% data_set)) {
    stop("data_set must be either K-array or H-array")
  }
  
  if("randomForest" %in% rownames(installed.packages()) == FALSE) {install.packages("randomForest")}
  library(randomForest)
  source("support_functions/load_AA_property_data.R") # loads the physical property data for the specified dataset
  
  # load data and response (signal to noise)
  temp_list = load_AA_property_data(data_set = data_set, position = c("X", "Z"), random = random)
  data_mat = temp_list[[1]]
  response = temp_list[[2]]
  
  # create model_func to interface with CV (cross validation function)
  model_func <- function(data_mat, response) {
    model_var = randomForest(data_mat, response, ntree = 1000, replace = T, localImp = F)
    return(model_var)
  }
  
  # make predict interface to return the intended test prediction by itself
  predict_func <- function(model, test_data) {
    test_prediction = predict(model, test_data, ntree = 1000, replace = T, localImp = F)
    return(test_prediction)
  }
  
  # random response interface.  Made for Q2 calculation
  random_response_func <- function(train_response, test_response) {
    random_response = rep(mean(train_response), length(test_response))
    return(random_response)
  }
  
  # performance function interface.  Can specify arbitrary performance function but default here is Q2
  source("ML_functions/Q2.R")
  performance_func <- function(train_response, test_prediction, test_response, performance_func_args = NULL) {
    
    continuous = !is.factor(unlist(test_response))
    # set variables used for make_random_response.  Can be specified with performance_func_args
    if (is.null(performance_func_args)) {
      if (continuous) {
        mean_random_prediction = mean(train_response)
        frac_correct = NULL
      } else {
        mean_random_prediction = NULL
        frac_correct = 1 / length(unique(train_response))
      }
    } else { # specified performance_func_args on the outside
      mean_random_prediction = performance_func_args$mean_random_prediction
      frac_correct = performance_func_args$frac_correct
    }
    
    # make random_response
    source("ML_functions/make_random_response.R")
    return_list = make_random_response(continuous = continuous, mean_random_prediction = mean_random_prediction, frac_correct = frac_correct)
    random_response_func = return_list$random_response_func
    random_response_func_args = return_list$random_response_func_args
    
    random_response = as.matrix(do.call(random_response_func, c(list(train_response, test_response), random_response_func_args)))
    
    return(Q2(test_prediction, test_response, random_response))
  }
  
  # model parameters
  k_folds = 5
  
  model_CV_obj = list(data_mat = data_mat, response = response, model_func = model_func, k_folds = k_folds, predict_func = predict_func, random_response_func = random_response_func, performance_func = performance_func)
  
  return(model_CV_obj)
  
}
