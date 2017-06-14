make_random_response <- function(continuous = T, mean_random_prediction = NULL, frac_correct = NULL) {
# makes standard random_response_func
# mean_random_prediction used for continuous.  NULL uses mean training response
# frac_correct used for classification

# make random_response_func
if (!continuous) {
  random_response_func <- function(train_response, test_response, frac_correct = NULL) {

    if (is.null(frac_correct)) { # this assumes the classes are relatively balanced
      frac_correct = 1 / length(unique(test_response))
    }
  
    test_response = as.factor(test_response)
  
    # scramble to get frac_correct accuracy
    random_response = test_response
    random_indices = sample(1:length(test_response), round(length(test_response) * (1 - frac_correct)))
  
    for (n_ind in random_indices) {
      other_label = levels(test_response)[!(levels(test_response) %in% test_response[[n_ind]])]
      random_response[[n_ind]] = other_label[[sample(length(other_label), 1)]]
    }
    return(as.matrix(random_response))
  }
  random_response_func_args = list(frac_correct = frac_correct)

} else { #######################################################

  random_response_func <- function(train_response, test_response, mean_random_prediction = NULL) {

    if (is.null(mean_random_prediction)) {
      mean_random_prediction = mean(train_response)
    }
    random_response = rep(mean_random_prediction, length(test_response))
    return(as.matrix(random_response))

  }
  random_response_func_args = list(mean_random_prediction = mean_random_prediction)

}

return(list(random_response_func = random_response_func, random_response_func_args = random_response_func_args))

}
