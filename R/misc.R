

# Normalize a matrix (m X n) by the sum of its elements
normalizeC <- function(x) x/.colSums(x, m = length(x), n = 1)  # a much faster version compared to sum or .colSums(x, m=1, n=length(x))

# Normalize a matrix (m X n) by rowSums
normalizeR <- function(x, m, n) x/.rowSums(x, m, n)  

order_states <- function(param) {
    ## ordering of estimates parameters per state, to force an ordering in the HMM
    ## param : n x Q matrix, where Q is the number of states return : n x Q where for
    ## each row a_i,q-1 < a_i,q , for q = 1, ..., Q
    aux = which(apply(param, 1L, function(r) any(diff(r) < 0L)))
    sorted = NULL
    if (length(aux)) 
        sorted = t(apply(param[aux, , drop=F], 1L, order))
    
    list(toSort = aux, order = sorted)
}
