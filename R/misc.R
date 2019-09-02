#' Miscellaneous functions for package traveltimeHMM
#' 
#' @name misc
NULL


#' Normalize a vector by the sum of its elements
#' @keywords internal
#' Function \code{normalizeV} normalizes a vector by the sum of its elements. It is much faster than sum or .colSums(x, m=1, n=length(x)) .
#' @param x  A vector of dimension \code{n}. 
#' @return A vector of dimension \code{n} with  its elements summing up to 1.
#' @examples
#' normalizeV(c(1,3))
#' 
#' \dontrun{
#' normalizeV(matrix(c(1,2,3,4), nrow=2, ncol = 2), 2, 2)
#' }
#' @export
normalizeV <- function(x) {
  if(class(x)=="numeric")
    x/.colSums(x, m = length(x), n = 1)  # a much faster version compared to sum or .colSums(x, m=1, n=length(x))
  else
    stop("Function normalizeV requires a vector.")
}

#' Normalize each row of a matrix by the sum of its elements
#' @keywords internal
#' Function \code{normalizeR} normalizes each row of a matrix by the sum of its elements.
#' @param x  A matrix of dimensions \code{m X n}.
#' @param m The number of rows of the matrix.
#' @param n The number of columns of the matrix.
#' @return A matrix of dimensions \code{m X n} with each row summing up to 1.
#' @examples
#' normalizeR(matrix(c(1,2,3,4), nrow=2, ncol = 2), 2, 2)
#' @export
normalizeR <- function(x, m, n) x/.rowSums(x, m, n)  


#' Order manually the states vector
#' @keywords internal
#' Function \code{order_states} orders each row of a states vector based on the mean.
#' @param param ...
#' @return A list of indeces to sort by and how to sort them
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
