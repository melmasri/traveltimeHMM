#' Computing the posterior marginals of hidden state variables using a sequence of observations
#' 
#' \code{forwardback} implements Russell and Norvig's forward-backward algorithm that is used
#' at step 3 of Algorithm 1 in Woodard et al.
#' 
#' @param probTran An \code{nObs x nQ^2} matrix having each row represent 
#' \eqn{P(Obs_k | state_k ) * P(state_k | state_{k-1})}.  Each row of \code{probTran} is such that
#' \code{matrix(probTran[1,], ncol=nQ, nrow=nQ, byrow =TRUE)} is the probability above in matrix format.
#' @param initQ A \code{nQ}-sized vector that represents the initial probability of the first link in a trip.
#' 
#' @return \code{probTran} returns a list of \code{alpha} and \code{beta}
#' (respectively forward and backward probability row-wise vectors) such
#' that \code{normalize(alpha[i] * beta[i])} is the smoothed HMM state estimate.
#' 
#' @references
#' {Russell, S., Norvig, P., 2009.  Artificial Intelligence: A Modern Approach.  Pearson Education, UK.}
#' 
#' @export


forwardback <- function(probTran, initQ) {
    fwd <- function(prob, i) normalizeC(.rowSums(probTran[i, ] * prob[indf], m = nQ, 
        n = nQ))
    bwd <- function(i, prob) normalizeC(.colSums(probTran[i, ] * prob[indb], m = nQ, 
        n = nQ))
    n = length(probTran)
    nQ = length(initQ)
    Obs = n/nQ^2
    indb = rep(1:nQ, nQ)
    indf = rep(1:nQ, each = nQ)
    if (!is.matrix(probTran)) 
        probTran = matrix(probTran, ncol = nQ^2, byrow = TRUE)
    ## probTranList = lapply(split(probTran, 1:Obs), function(r) matrix(r, ncol=nQ,
    ## nrow=nQ, byrow=TRUE))
    list(alpha = unlist(Reduce(fwd, 1:Obs, init = initQ, acc = TRUE)[-1], use.names = FALSE), 
        beta = unlist(Reduce(bwd, 1:Obs, init = matrix(rep(1, nQ), nrow = nQ), acc = TRUE, 
            right = TRUE)[-1], use.names = FALSE))
}

