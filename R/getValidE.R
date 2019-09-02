#' Get an appropriate value for \code{logE}
#' @keywords internal
#' 
#' Function \code{getValidE} takes a tentative value for \code{logE} and returns a valid value for it.
#' 
#' The value for \code{logE} in the model object passed as parameter is \emph{NOT} used.
#' 
#' @param object A model object (a list) provided through the execution of function \code{timetravelHMM}.
#' @param logE Tentative values of point estimates of trip effects.
#' @param n Number of samples. Default is 1000.
#' 
#' @return The function return a valid value for \code{logE}.  The value supplied as a parameter
#'   for \code{logE} is a numerical vector of size \code{n}, in which case it is
#'   returned \emph{as is}.  If a single numeric value is supplied, it will be replicated into
#'   a vector of the appropriate size which will be returned.  If \code{logE} is \code{NULL}
#'   or invalid, then the function will return either a vector of simulated values
#'   (if the model is from the \code{trip} family), or a vector of \code{0} otherwise.
#' @export
getValidE <- function(object, logE, n) {
  # Part 1 of validation.  We check for the existence of a valid parameter logE
  # of size either 1 or n
  if(!is.null(logE)) { # If a numeric parameter exists only...
    
    # If parameter is a single numeric, convert to vector of size n
    if(is.numeric(logE) && length(logE) == 1)
      logE <- rep(logE, n)
    
    # If we don't end up with a vector of size n, send warning message.
    if(!is.numeric(logE) || length(logE) != n) {
      message("Values in parameter logE need to be a vector of length 'n'.  Values for logE are discarded.")
      logE <- NULL
    }
  }
  # At this point, logE has the appropriate value among valid parameters passed,
  # and NULL if no such value is found.
  
  # Part 2 of validation: if logE is null...
  
  if(is.null(logE)) {
    if(grepl('trip', object$model)) 
      logE = rnorm(n, mean = 0, sd = object$tau) # ... we define a vector of random values if trip model...
    else logE = rep(0, n) # ... otherwise we set to vector of zeros.
  }
  logE # We return the result.
  # How do we handle tau if logE is set to 0?  Currently it keeps its value.  EG 2019/07/25
}
