#' Commputes initial state probabilites 
#'
#' \code{initial_est} return the state probabilities for each level in the passed \code{by_factor},
#' defined as the expectation of the vectors for each level.
#'
#' @param state_prob An \code{n x Q} matrix of initial state probabilities for \code{n} observations and \code{Q} states
#' @param by_factor A vector of factors of size \code{n} (a factor for each row of \code{state_prob}).
#' @param subset A vector indexing the rows of \code{state_prob} for which states probabilities are computed using only this subset, defualt is \code{NULL}
#'
#' @details NULL
#' 
#' @return An \code{m x Q} matrix of probabilites for the m levels of \code{by_factor}.
#'
#' 
#' @examples
#' x = runif(10)
#' x = cbind(x, 1-x)
#' initial_est(x, factor(sample(c("A", "B"), 10, replace=TRUE)))
initial_est <- function(state_prob, by_factor, subset=NULL) {
    # Basic verfication of inputs; stop if any is incorrect.
    if(!is.matrix(state_prob))
        stop('state_prob is not a matrix!')
    if(!is.factor(by_factor))
        stop('by_factor is not a factor')
    if(nrow(state_prob) != length(by_factor))
        stop('legnth of factors is not the same is the number of rows of state_prob.')
    if(!is.null(subset) & !is.vector(subset)) # If not null, 'subset' needs to be a vector.
        stop('subset must be  vector.')

    nQ = ncol(state_prob) # Number of states as defined elsewhere
    
    # Compute a matrix of dimensions nQ x m, either using the full set of observations
    # or the subset defined in "subset". 
    if(is.null(subset)){
        init = vapply(split(state_prob, by_factor),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE), numeric(nQ), USE.NAMES = FALSE)
    }else{
        init = vapply(split(state_prob[subset, ], by_factor[subset]),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE), numeric(nQ), USE.NAMES = FALSE)
    }
    
    # return the final, starightforward nQ x m matrix
    return(matrix(drop(init), ncol = nQ, byrow = TRUE))
}


#' @describeIn Commputes transition matrix probabilities
#'
#' \code{tmat_est} return transition state probabilities for each level in the passed \by{by_factor}.
#'
#' @param joint_prob An \code{n x Q^2} matrix of joint state probabilities for \code{n} observations and \code{Q} states
#' @param state_prob An \code{n x Q} matrix of initial state probabilities for \code{n} observations and \code{Q} states
#' @param init_ids A vector of indices of the initial state observation for each trip
#' @param by_factor A vector of factors of size \code{n} (a factor for each row of \code{state_prob}).

#'
#' @return An \code{m x Q} matrix of probabilites for the m levels of \code{by_factor}.
#'
#' @examples
#' \dontrun{
#' }
tmat_est <- function(joint_prob,state_prob, init_ids, by_factor){
    # Basic verfication of inputs; stop if any is incorrect.
    if(!is.matrix(joint_prob))
      stop('state_prob is not a matrix!')
    if(!is.matrix(state_prob))
      stop('state_prob is not a matrix!')
    if(!is.factor(by_factor))
      stop('by_factor is not a factor')
    if(nrow(joint_prob) != length(by_factor)-1)
      stop('legnth of factors does not match the number of rows of joint_prob.')
    if(nrow(state_prob) != length(by_factor))
      stop('legnth of factors is not the same is the number of rows of state_prob.')
    nQ <- ncol(state_prob) # Number of states
    nQ2 <- nQ^2 # Compute nQ2, the square of nQ, which will be used often
    nObs <- nrow(state_prob) # Number of observations

    # We need to compute a ratio.  (The first observation is removed from the factors.)
    # Numerator (num): compute a matrix of dimensions nQ^2 x m (where m = the number of factors)...
    num = vapply(split(joint_prob, by_factor[-1]), function(r) .colSums(r,length(r)/nQ2, nQ2), numeric(nQ2), USE.NAMES = FALSE)
    num = matrix(drop(num), ncol = nQ2, byrow = TRUE) # ...which becomes a cleaner nQ^2 x m matrix.
    
    # Denominator (den)
    den <- state_prob # Start with the initial state probability matrix...
    den[init_ids[-1] - 1, ] <- 0 # Setting the overlapping multiplication between trips to 0
    
    # Compute a matrix of dimensions nQ x m (where m = the number of factors)...
    den = vapply(split(den[-nObs, ], by_factor[-1]), function(r) .colSums(r, length(r)/nQ, nQ), numeric(nQ), USE.NAMES = FALSE)
    den = matrix(drop(den), ncol = nQ, byrow = TRUE) # ...which becomes a cleaner nQ x m matrix.
    
    # We return the ratio of 'num' over a modified 'den' with each column replicated nQ times consecutively.
    return(num/den[, rep(1:nQ, each = nQ)])
}
