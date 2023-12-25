
# -------------------------------------------------------------------------
#
# function to analyze estimated ECV
#
# -------------------------------------------------------------------------

## Arguments

# object: SimResult object being described.
# improper: If TRUE, include the replications that provided improper solutions.
# digits: The number of digits rounded in the result. If NULL, the results will not be rounded.


## Value

# Population Value: population value of the ECV
# Estimate Average: average estimated ECV across all replications
# Estimate SD: standard deviation of estimated ECV across all replications
# Average Bias: Difference between estimated ECV and population ECV
# Relative Bias: relative bias, which is (Estimate Average - Population Value)/Population Value



summaryECV <- function(object, improper = TRUE, digits = NULL){

  # population parameters ---------------------------------------------------
  
  # std parameters
  popParams <- object@stdParamValue
  
  # only loadings
  popLoadings <- popParams[, grep('=~', names(popParams))]
  
  # determine g factor
  facLabels <- sub('=~.*', "", names(popLoadings))
  facLabelG <- names(which.max(table(facLabels))) # g = factor with largest number of indicators
  
  # ECV in population
  popLoadingsG <- popLoadings[, grepl(facLabelG, names(popLoadings))]
  popLoadingsS <- popLoadings[, !grepl(facLabelG, names(popLoadings))]
  popECV <- sum(popLoadingsG^2) / sum(popLoadingsG^2, popLoadingsS^2)
  
  
  # estimated parameters ----------------------------------------------------
  
  # std estimates
  if(improper == TRUE){
    estParams <- object@stdCoef
  } else{
    estParams <- object@stdCoef[object@converged == 0, ]
  }
  
  # only loadings
  estLoadings <- estParams[, grep('=~', names(estParams))]
  
  # estimated ECV
  getECV <- function(x){
    num <- sum(x[grepl(facLabelG, names(x))]^2)
    denom <- sum(x^2)
    ECV <- num/denom
    return(ECV)
  }
  
  estECV <- apply(estLoadings, 1, getECV)
  
  
  # results -----------------------------------------------------------------
  
  average.est <- mean(estECV, na.rm = T)
  sd.est <- sd(estECV, na.rm = T)
  average.bias <- mean(estECV - popECV, na.rm = T)
  relative.bias <- (average.est - popECV)/popECV
  
  result <- cbind(popECV, average.est, sd.est, 
                  average.bias, 
                  relative.bias)
  
  colnames(result) <- c("Population Value", "Estimate Average", "Estimate SD", 
                        "Average Bias",
                        "Relative Bias")
  rownames(result) <- 'ECV'
  
  if(is.null(digits)){
    result
  } else{
    round(result, digits = digits)
  }

}








