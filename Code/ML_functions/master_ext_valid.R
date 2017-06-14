master_ext_valid <- function(ext_valid_obj, return_prediction = F) {
  # this function performs external cross validation
  
  source("ML_functions/Q2.R")
  
  #### unpack the variables
  train_data = ext_valid_obj$train_data
  test_data = ext_valid_obj$test_data
  train_response = ext_valid_obj$train_response
  test_response = ext_valid_obj$test_response
  model_func = ext_valid_obj$model_func
  model_func_args = ext_valid_obj$model_func_args
  predict_func = ext_valid_obj$predict_func
  predict_func_args = ext_valid_obj$predict_func_args
  random_response_func = ext_valid_obj$random_response_func
  random_response_func_args = ext_valid_obj$random_response_func_args
  
  # use standard predict function if none was specified
  if (is.null(predict_func)) {
    predict_func = predict
  }
  
  ##### training and testing phase
  model = do.call(model_func, c(list(train_data, train_response), model_func_args))
  test_prediction = do.call(predict_func, c(list(model, test_data), predict_func_args))
  random_response = do.call(random_response_func, c(list(train_response, test_response), random_response_func_args))
  
  # wrap in list
  original_names = names(test_prediction)
  test_prediction = lapply(1:length(test_prediction), function(i) test_prediction[[i]])
  names(test_prediction) = original_names
  
  # calculate Q2
  Q2 = Q2(test_prediction, test_response, random_response)
  
  # rare errors during model testing creates nan Q2
  if (is.nan(Q2) || is.na(Q2)) {
    Q2 = 0
  }
  
  if (return_prediction) {
    return(list(Q2, test_prediction))
  } else {
    return(Q2)
  }
  
}
