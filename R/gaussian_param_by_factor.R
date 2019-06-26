#' @keywords internal
#' Calculate the parameters of Gaussian speeds
#' 
#' \code{gaussian_param_by_factor} calculates the mean and variance of Gaussian speeds.  The function refers to the calculation
#' of \eqn{\mu} and \eqn{\sigma^{2}} at step 4 of Algorithm 1 in Woodard et al., 2017
#' 
#' @param logspeeds A numeric vector of speed observations (in km/h) on the (natural) log-scale, of size \code{n = number of observations}.
#' @param byFactor A vector of factors (link + time bin) of size \code{n}, i.e. one factor for each speed observation.
#' @param state_prob An \code{n x nQ} matrix of smoothed state probabilities where \code{nQ = number of states}.
#' @details NULL
#' 
#' @return \code{gaussian_param_by_factor} returns a list of two matrices of identical dimensions \code{m X nQ}
#' where \code{m = number of levels} (links + timeBins).  The first list \code{mean} contains the vector of means
#' for each level and state.  The second list \code{sigma2} contains the vector of variances for each level and state.  
#' @references
#' {Woodard, D., Nogin, G., Koch, P., Racz, D., Goldszmidt, M., Horvitz, E., 2017.  Predicting travel time reliability using mobile phone GPS data.  Transportation Research Part C, 75, 30-44.}
#' @export

gaussian_param_by_factor <- function(logspeeds, byFactor, states = NULL) {
    if (is.null(states)) { # If 'states' is not provided, we assume a single state '1' and proceed.
        nQ = 1
        states = 1
        sumStates = unlist(lapply(split(rep(1, length(logspeeds)), byFactor), function(r) .colSums(r, 
            length(r)/nQ, nQ)), use.names = FALSE)
    } else {
        nQ = ncol(states)
        sumStates = unlist(lapply(split(states, byFactor), function(r) .colSums(r, 
            length(r)/nQ, nQ)), use.names = FALSE) # We get a m X nQ matrix of the sum of phi_j,k
                                                   # for each level, corresponding to the denominator
                                                   # of both equations.
    }
    sumSpeed = unlist(lapply(split(states * logspeeds, byFactor), function(r) .colSums(r, 
        length(r)/nQ, nQ)), use.names = FALSE) # In the same way we compute the numerator
                                               # of the mu equation (sum of (states * logspeeds)...
    sumSpeed2 = unlist(lapply(split(states * (logspeeds^2), byFactor), function(r) .colSums(r, 
        length(r)/nQ, nQ)), use.names = FALSE) # ... and part of that of the sigma^2 equation
                                               # (sumSpeed2 = sum of (states * logspeeds^2)
    
    mean = sumSpeed/sumStates # compute mu : divide numerator by sumStates
    sig2 = sumSpeed2/sumStates - mean^2 # compute sigma^2 as E(logspeeds^2) - mu^2
                                        # where E(logspeeds^2) = sumSpeed2 / sumStates
                    
    list(mean = matrix(mean, ncol = nQ, byrow = TRUE),
         sigma2 = matrix(pmax(sig2, 1e-10), ncol = nQ, byrow = TRUE))
}
