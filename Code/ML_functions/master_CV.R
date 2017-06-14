master_CV <- function(model_CV_obj) {
# this function does cross validation on given train data and response.  Must specify a model_func to interface with the CV function
# model_func must take the form of model_func(data_mat, train_response)
# this function returns a CV_results list

source("ML_functions/CV.R")

# unpack the variables
data_mat = model_CV_obj$data_mat
response = model_CV_obj$response
model_func = model_CV_obj$model_func
k_folds = model_CV_obj$k_folds
predict_func = model_CV_obj$predict_func
performance_func = model_CV_obj$performance_func

# if custom arguments are specified, set them.  Else they are NULL
model_func_args = NULL
if ("model_func_args" %in% names(model_CV_obj)) {
  model_func_args = model_CV_obj$model_func_args
}
predict_func_args = NULL
if ("predict_func_args" %in% names(model_CV_obj)) {
  predict_func_args = model_CV_obj$predict_func_args
}
performance_func_args = list()
if ("performance_func_args" %in% names(model_CV_obj)) {
  performance_func_args = model_CV_obj$performance_func_args
}

model = do.call(model_func, c(list(data_mat), list(response), model_func_args))

##### evaluate the model with current feature set #####
return_list = CV(data_mat = data_mat, response = response, k_folds = k_folds, 
  model_func = model_func, model_func_args = model_func_args, 
  predict_func = predict_func, predict_func_args = predict_func_args, 
  performance_func = performance_func, performance_func_args = performance_func_args)

return(c(return_list, model = list(model)))

}


