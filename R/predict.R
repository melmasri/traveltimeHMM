
#' @export
predict.traveltime<-function(object, data, starttime = Sys.time(),  n=1000, ...){
    if(!is.list(data))
        stop('data must be a list, data.fram or data.table')
    if(!all(c('linkId', 'length') %in% names(data)))
        stop('data must have objects named linkId and length, corresponding order of travelled links and their lengths')
    if(length(data$linkId)!=length(data$length))
        stop('length of objects do not match!')
    
    if(object$model == 'no-dependence')
        return(predict.traveltime.no_dependence(object, data , starttime, n, ...))
    
    if(object$model =='trip')
        return(predict.traveltime.no_dependence(object, data , starttime, n, ...))
    
    if(grepl('HMM', object$model))
        return(predict.traveltime.HMM(object, data, starttime, n, ...))
}

predict.traveltime.no_dependence <- function(object, data, starttime, n = 1000, ...) {
    param = list(...)
    linkIds = data$linkId
    len = data$length
    ## sampling E (trip-effect)
    if(!is.null(param$E) && is.numeric(param$E)){
        E = param$E
    }else if(grepl('trip', object$model))
        E = rnorm(n, mean = 0, sd = object$tau) else E = 0

    fact = paste(linkIds[1], time_bins(starttime), sep = ".")
    id = which(levels(object$factors) == fact)
    speed = rnorm(n, object$mean[id, ], object$sd[id, ])
    tt = len[1] * exp(-speed - E)
    for (k in 2:length(linkIds)) {
        fact = as.factor(paste(linkIds[k], time_bins(starttime + tt), sep = "."))
        id = sapply(levels(fact), function(s) which(levels(object$factors) == s, useNames = FALSE), USE.NAMES = FALSE)
        ind = as.numeric(fact)
        speed = rnorm(n, object$mean[id[ind], ], object$sd[id[ind], ])
        tt = tt + len[k] * exp(-speed - E)
    }
    tt
}


predict.traveltime.HMM <- function(object, data, starttime, n = 1000, ...) {
    ## sampling E (trip-effect)
    param = list(...)
    linkIds = data$linkId
    len = data$length
    if(!is.null(param$E) && is.numeric(param$E)){
        E = param$E
    }else if(grepl('trip', object$model))
        E = rnorm(n, mean = 0, sd = object$tau) else E = 0
    nQ = object$nQ
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
                         sep = "."), useNames = FALSE)
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    tt = len[1] * exp(-speed - E)
    if (length(linkIds) > 1) 
        for (k in 2:length(linkIds)) {
            ## this saves about 50ms for 117 routes (.4 ms per route)
            fact = as.factor(paste(linkIds[k], time_bins(starttime + tt), sep = "."))
            id = sapply(levels(fact), function(s) which(levels(object$factors) == 
                s, useNames = FALSE), USE.NAMES = FALSE)
            ind = as.numeric(fact)
            indIds = (1:length(id) - 1) * nQ
            ## indIds is suppsed to be the first index of each unique factor, if length(id) =1 , indIds = c(0,0), if length(id)=2 then inIds = c(0,2)
            ## the logic behind indIds, is that if nQ =2, then the first two rows for tmat2 corresponds to the first unique factor, the second 2 rows to the second unique factor
            ## etc, .. the number of unique factors is length(id)
            ## creating tmat takes about 10ms for 2 rows or the full matrix to speed up create
            ## tmat2 before the loop or post-estimation for the whole Quebec dataset, for 1294
            ## obs it takes about an extra 11sec if at the top of the loop
            tmat2 = matrix(c(t(object$tmat[id, ])), ncol = nQ, byrow = TRUE)
            ## tmat2 is an nQ x length(id)
            ## the first two rows are the transition of levels(fact)[1], the seond 2 rows are the transitions of levels(fact)[2], and so on.
            # if nQ =2 , lenth(id)=1 this is 2x2 matrix tmat[1, ] = (q1,q1), (q1,q2), (q2,q1), (q2,q2),
            ## tmat2  = (q1,q1) , (q2, q1), 2nd row (q1,q2), (q2,q2), tmat2 has rows equal to unique factors of fact
            tmat2 = t(apply(tmat2, 1, cumsum))                                 # here tmat2 = c(total prob. to go to q1, total prob. to go to q2) = c(p(q1), p(q2)) if length(id)=1 an nQ=2
            tmat2 = tmat2[indIds[ind] + Qk, ]                                  # QK = 1, 2, 1, 2, 2, 1, 2, ... is the old states
            ## since indIds is the index of the top row of the transition matrix for each factor
            ## then, asumme ind correspond to the 2nd unique factor
            ## then inIds[ind] would return 2, meaning row 3, and 4 of tmat2 is transition matrix of this fact
            ## if the old state Qk=1, then 2+1, returns row 3, if old state is 2, this returns row 4.
            ## Another example, assume that ind = 4th unique factor, then indIds[4] = 3*2 = 6, meaning the 4th factor has transition matrix in levels 7 and 8.
            ## if the old state is 1 then  tmat2[indIds[ind] + Qk, ] returne row 7 = indIds[4] + 1, otherwise returns row 8
            ## and so on.
            ## QK below samples a new state based on the new tmat2, which correspods now to 1 row per unqiue factor, and it is the row of the transition matrix of the old state.
            ## now that I think about it, I am not sure about this line tmat2 = t(apply(tmat2, 1, cumsum)), need to be checked.
            Qk = max.col(runif(n) < tmat2, "first")
            speed = rnorm(n, object$mean[cbind(id[ind], Qk)], object$sd[cbind(id[ind], Qk)])
            tt = tt + len[k] * exp(-speed - E)
        }
    tt
}

