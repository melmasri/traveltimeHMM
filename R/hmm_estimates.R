#' Commputes initial state probabilites 
#'
#' \code{initial_est} return the state probabilities for each level in the passed \code{by_factor}.
#'
#' @param state_prob An \code{n x Q} martrix of initial state probabilities for \code{n} observations and \code{Q} states
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
    if(!is.matrix(state_prob))
        stop('state_prob is not a matrix!')
    if(!is.factor(by_factor))
        stop('by_factor is not a factor')
    if(nrow(state_prob) != length(by_factor))
        stop('legnth of factors is not the same is the number of rows of state_prob.')
    if(!is.null(subset) & !is.vector(subset))
        stop('subset must be  vector.')

    nQ = ncol(state_prob)
    if(is.null(subset)){
        init = vapply(split(state_prob, by_factor),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE), numeric(nQ), USE.NAMES = FALSE)
    }else{
        init = vapply(split(state_prob[subset, ], by_factor[subset]),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE), numeric(nQ), USE.NAMES = FALSE)
    }
    return(matrix(drop(init), ncol = nQ, byrow = TRUE))
}


#' @describeIn Commputes transition matrix probabilities
#'
#' \code{tmat_est} return transition state probabilities for each level in the passed \by{by_factor}.
#'
#' @param joint_prob An \code{n x Q} martrix of initial state probabilities for \code{n} observations and \code{Q} states
#' @param state_prob An \code{n x Q} martrix of initial state probabilities for \code{n} observations and \code{Q} states
#' @param init_ids 
#' @param by_factor A vector of factors of size \code{n} (a factor for each row of \code{state_prob}).

#'
#' @return An \code{m x Q} matrix of probabilites for the m levels of \code{by_factor}.
#'
#' @examples
#' \dontrun{
#' }
tmat_est <- function(joint_prob,state_prob, init_ids, by_factor){
    nQ <- ncol(state_prob)
    nQ2 <- nQ^2
    nObs <- nrow(state_prob)

    joint_prob = joint_prob[-c(init_ids[-1] - 1), ] # removing first obs per trip
    by_factor  = by_factor[-init_ids]               # removing first obs per trip
    
    num = vapply(split(joint_prob, by_factor), function(r) .colSums(r,length(r)/nQ2, nQ2), numeric(nQ2), USE.NAMES = FALSE)
    num = matrix(drop(num), ncol = nQ2, byrow = TRUE)
    
    den <- state_prob[-c(init_ids[-1] - 1, nObs),] # removing last obs per trip
    den = vapply(split(den, by_factor), function(r) .colSums(r, length(r)/nQ, nQ), numeric(nQ), USE.NAMES = FALSE)
    den = matrix(drop(den), ncol = nQ, byrow = TRUE)
    return(num/den[, rep(1:nQ, each = nQ)])
}

