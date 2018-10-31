
forwardback<-function(probTran, initQ){
    ## initQ is the initial probability of the first link in a trip, it is a nQ-vector
    ## probTran is an nObs x nQ^2 matrix with:
    ## each row of probTran reprsents P(Obs_k | state_k ) * P(state_k | state_{k-1})
    ## ProbTran is a row-wise matrix such that matrix(probTran[1,], ncol=nQ, nrow=nQ, byrow =TRUE) is the probability above in matrix format
    ## returns a list of alpha = forward prob,  beta = backward prob, such that normalize(alpha[i] *beta[i]) is the smoothed HMM state estimate
    ## alpha and beta are row-wise matrix vectors
    fwd <-function(prob, i)  normalizeC(.rowSums(probTran[i,] * prob[indf], m=nQ, n=nQ))
    bwd <-function(i, prob)  normalizeC(.colSums(probTran[i,] * prob[indb], m=nQ, n=nQ))
    n = length(probTran)
    nQ = length(initQ)
    Obs = n/nQ^2
    indb = rep(1:nQ, nQ)
    indf = rep(1:nQ, each = nQ)
    if(!is.matrix(probTran)) probTran = matrix(probTran, ncol=nQ^2, byrow=TRUE)
    ## probTranList = lapply(split(probTran, 1:Obs), function(r) matrix(r, ncol=nQ, nrow=nQ, byrow=TRUE))
    list(alpha = unlist(Reduce(fwd, 1:Obs, init = initQ, acc = TRUE)[-1], use.names=FALSE),
         beta  = unlist(Reduce(bwd, 1:Obs, init = matrix(rep(1, nQ), nrow=nQ), acc = TRUE, right = TRUE)[-1],use.names=FALSE))
}

