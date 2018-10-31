
#' Commputes initial state probabilites 
#' 
#' @param state_prob An n x Q martrix of initial 
#' @param by_factor A vector of factors of size n (a factor for each row of \code{state_prob}).
#' @param sunset A vector 
#' @return 
#' @examples
#' \dontrun{add(1, 1)}
#' 

initial_est <- function(state_prob, by_factor, subset=NULL, second_factor=NULL) {
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
        init = lapply(split(state_prob, by_factor),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE))
    }else{
        init = lapply(split(state_prob[subset, ], by_factor[subset]),
            function(r) .colMeans(r,length(r)/nQ, nQ, na.rm = TRUE))
    }
    init = matrix(unlist(init), ncol = nQ, byrow = TRUE)

    if(!is.null(second_factor)){
        init_impute = initial_est(state_prob, second_factor, subset)
        
    }
    
        
}


