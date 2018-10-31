

gaussian_param_by_factor <- function(speeds, byFactor, states = NULL) {
    ## calculating mean and variance of gaussian speeds : vector of logspeeds states :
    ## an n X nQ matrix of smoothed state probabilities. nQ No. of states, n No. of
    ## Obs byFactor : a vector of factors of the same length as speeds, with with a
    ## factor for each speed observation returns a list of meean snd sigma of Gaussian
    ## parameters per Obs.
    if (is.null(states)) {
        nQ = 1
        states = 1
        sumStates = unlist(lapply(split(rep(1, length(speeds)), byFactor), function(r) .colSums(r, 
            length(r)/nQ, nQ)), use.names = FALSE)
    } else {
        nQ = ncol(states)
        sumStates = unlist(lapply(split(states, byFactor), function(r) .colSums(r, 
            length(r)/nQ, nQ)), use.names = FALSE)
    }
    sumSpeed = unlist(lapply(split(states * speeds, byFactor), function(r) .colSums(r, 
        length(r)/nQ, nQ)), use.names = FALSE)
    sumSpeed2 = unlist(lapply(split(states * (speeds^2), byFactor), function(r) .colSums(r, 
        length(r)/nQ, nQ)), use.names = FALSE)
    
    mean = sumSpeed/sumStates
    sig = sumSpeed2/sumStates - mean^2
    list(mean = matrix(mean, ncol = nQ, byrow = TRUE), sigma = matrix(pmax(sig, 0), 
        ncol = nQ, byrow = TRUE))
}
