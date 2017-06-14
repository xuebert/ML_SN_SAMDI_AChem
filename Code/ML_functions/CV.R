CV <- function(data_mat, response, model_func, model_func_args = NULL, k_folds = nrow(data_mat), predict_func = NULL, predict_func_args = NULL, performance_func, performance_func_args = NULL, ...) {
  # this function automatically does cross-validation on model_func and saves the results as a list 
  # this function is for internal cross validation and should only be called by master_CV
  
  num_rows = nrow(data_mat)
  
  # LOO is a special case that does not do random k-folds.  Just keep the indices in order
  if (k_folds == num_rows) {
    test_indices = 1:num_rows
  } else {
    
    # function that divides a vector into n equal sized chunks
    chunk <- function(x, n) {split(x, factor(sort(rank(x) %% n)))}
    
    # make equal sized chunks for k-folds
    test_indices = chunk(sample(1:num_rows, num_rows), k_folds)
    
  }
  
  # storage variables
  test_prediction = vector("list", length(test_indices))
  test_response = test_prediction
  counter = 0
  
  if (is.null(predict_func)) {
    predict_func = predict
  }
  
  # convert response to a matrix
  if (is.vector(response)) {
    response = matrix(response, ncol = 1, dimnames = list(names(response)))
  }
  
  performance = c()
  for (n_test in test_indices) {
    counter = counter + 1
    
    # specify datasets
    train_data = data_mat[-n_test, , drop = F]
    colnames(train_data) = colnames(data_mat)
    train_response = response[-n_test, , drop = F]
    test_data = data_mat[n_test, , drop = F]
    colnames(test_data) = colnames(data_mat)
    rownames(test_data) = rownames(data_mat)[n_test]
    
    # create model on training data
    model = do.call(model_func, c(list(train_data, train_response), model_func_args))
    
    # predict held out data set
    test_prediction[[counter]] = do.call(predict_func, c(list(model), list(test_data), predict_func_args))
    
    test_response[[counter]] = as.matrix(response[n_test,])
    
    performance = c(performance, do.call(performance_func, c(list(train_response, test_prediction[[counter]], test_response[[counter]]), performance_func_args)))
  }
  
  return(list(performance = performance, test_indices = test_indices, test_prediction = test_prediction, test_response = test_response))
  
}
