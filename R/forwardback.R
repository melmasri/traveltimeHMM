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
  
    # 'fwd' corresponds to the FORWARD function as defined in Eq. 15.5 of Russell and Norvig's book.
    # Argument 'prob' corresponds to fv[i-1] of figure 15.4 in the book, i.e. a single forward message
    # vector of dimension nQ for the steps up to the previous one. Argument 'i' is an integer which
    # represents an observation ID. Output corresponds to 'fv[i]', i.e. a single forward message vector
    # of dimension nQ.
    # 
    # WARNING: the function depends on the values of probTran, indf and nQ which are defined
    # one level above and ARE NOT passed as parameters.  This is suboptimal but necessary due
    # to the fact that fwd is called by the Reduce function which has strict requirements
    # on the function's interface.
  
    fwd <- function(prob, i) {
      return(normalizeV(.rowSums(probTran[i, ] * prob[indf], m = nQ, 
        n = nQ))) # Normalize by column.  This implementation is quite different from that in the paper,
                  # where some interaction happens between fv[i] and b.  ÉG 2019/06/10.
    }

    # 'bwd' corresponds to the BACKWARD function as defined in Eq. 15.9 of Russell and Norvig.  Arguments are reversed
    # (relative to 'fwd' and also to the paper) because of the call by Reduce with option 'right = TRUE'.
    # Argument 'prob' corresponds to b of figure 15.4 in the book, i.e. a single backward message
    # vector of dimension nQ for the steps down to the next one. Argument 'i' is an integer which
    # represents an observation ID. Output corresponds to 'b' a single backward message vector
    # of dimension nQ.
    # 
    # WARNING: the function depends on the values of probTran, indb and nQ which are defined
    # one level above and ARE NOT passed as parameters.  This is suboptimal but necessary due
    # to the fact that fwd is called by the Reduce function which has strict requirements
    # on the function's interface.
    #
    # QUESTION: 'bwd' does NOT make use of 'fv' from 'fwd' unlike in the paper.  The interaction between 'fv' and 'b'
    # seem to occur outside of function forwardback.
    
    bwd <- function(i, prob) {
      return(normalizeV(.colSums(probTran[i, ] * prob[indb], m = nQ, 
        n = nQ))) # TO DO: check that it does what is intended.  WHY do we want to normalize on BWD?
    }
    
    n = length(probTran)
    nQ = length(initQ) # Corresponds to the number of states
    Obs = n/nQ^2 # probTrans has dimensions [n/nQ^2, nQ^2], so Obs is the number of its rows.
    
    # 'indb' and 'indf' are used in the workings of functions 'fwd' and 'bwd'.
    indb = rep(1:nQ, nQ)
    indf = rep(1:nQ, each = nQ)

    # We make sure that probTran has the form of a matrix, and adjust accordingly.
    if (!is.matrix(probTran)) 
        probTran = matrix(probTran, ncol = nQ^2, byrow = TRUE)

    return(list(alpha = unlist(Reduce(fwd, 1:Obs, init = initQ, acc = TRUE)[-1], use.names = FALSE), 
        beta = unlist(Reduce(bwd, 1:Obs, init = matrix(rep(1, nQ), nrow = nQ), acc = TRUE, 
            right = TRUE)[-1], use.names = FALSE)))
}

