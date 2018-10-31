posterior_initial_states<-function(states, byFactor, index){
    ## computes colMeans of states, based on byFactor, and the subset indicated by index
    ## states      : n x nQ matrix, for Q states
    ## byFactor    : a factor vector of length n
    ## index       : an index of subset of 1, ..,n
    nQ = ncol(states)
    meanStates = lapply(split(states[index,], byFactor[index]), function(r) .colMeans(r, length(r)/nQ, nQ, na.rm=TRUE))
    a = matrix(unlist(meanStates), ncol=nQ, byrow=TRUE)
    rownames(a)<-levels(byFactor)
    a
}
