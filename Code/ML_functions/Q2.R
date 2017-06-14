Q2 <- function(test_prediction, test_response, random_response) {
# this function calculates Q2 values

# detect if list.  Assume each element in list if its own validation set.  if test_prediction is list, then assume all others are lists
if (class(test_prediction) == "list") {
  num_loop = length(test_prediction)
} else {
  num_loop = 1
  test_prediction = list(test_prediction)
  test_response = list(test_response)
  random_response = list(random_response)
}

# set up squared total error and squared residual error
SS_tot = 0
SS_res = 0

# iterate over each test_prediction
for (n_loop in 1:num_loop) {

  test_prediction_temp = test_prediction[[n_loop]]
  test_response_temp = test_response[[n_loop]]
  random_response_temp = random_response[[n_loop]]
  
  # continuous variable
  if (!is.factor(test_prediction_temp)) {
  
    # accumulate errors
    SS_tot = SS_tot + sum((test_response_temp - random_response_temp) ^ 2)
    SS_res = SS_res + sum((test_response_temp - test_prediction_temp) ^ 2)
    
  } else { # classification
  
    # predicted classes don't include all classes in test response
    if (!identical(levels(test_prediction_temp), levels(test_response_temp))) {
      levels(test_prediction_temp) = levels(test_response_temp)
    }

    # accumulate errors
    SS_tot = SS_tot + sum(test_response_temp != random_response_temp)
    SS_res = SS_res + sum(test_response_temp != test_prediction_temp)
  
  }
}

# do final Q2 calculation
Q2 = 1 - SS_res / SS_tot

# something went wrong and there was NA (rare)
if (is.nan(Q2) || is.na(Q2)) {
  Q2 = 0
  warning("Q2 in Q2.R created nan")
}

return(Q2)

}
