
predict.traveltime<-function(object, linkIds, len, starttime,  n=1000){
    if(object$model == 'no-dependence')
        return(predict.traveltime.no_dependence(object, linkIds, len , starttime, n))
    
    if(object$model =='trip')
        return(predict.traveltime.no_dependence(object, linkIds, len , starttime, n))
    
    if(grepl('HMM', object$model))
        return(predict.traveltime.HMM(object, linkIds, len, starttime, n))
}

predict.traveltime.no_dependence <- function(object, linkIds, len, starttime, n = 1000) {
     ## sampling E (trip-effect)
    if(grepl('trip', object$model))
        E = rnorm(n, mean = 0, sd = est$tau) else E = 0
    
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


predict.traveltime.HMM <- function(object, linkIds, len, starttime, n = 1000, ...) {
     ## sampling E (trip-effect)
    if(grepl('trip', object$model))
        E = rnorm(n, mean = 0, sd = object$tau) else E = 0
    nQ = object$nQ
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
        sep = "."), useNames = FALSE)
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    tt = len[1] * exp(-speed - E)
    if (length(linkIds) > 1) 
        for (k in 2:length(linkIds)) {
            # this saves about 50ms for 117 routes (.4 ms per route)
            fact = as.factor(paste(linkIds[k], time_bins(starttime + tt), sep = "."))
            id = sapply(levels(fact), function(s) which(levels(object$factors) == 
                s, useNames = FALSE), USE.NAMES = FALSE)
            ind = as.numeric(fact)
            indIds = (1:length(id) - 1) * nQ
            ## creating tmat takes about 10ms for 2 rows or the full matrix to speed up create
            ## tmat2 before the loop or post-estimation for the whole Quebec dataset, for 1294
            ## obs it takes about an extra 11sec if at the top of the loop
            tmat2 = matrix(drop(t(object$tmat[id, ])), ncol = nQ, byrow = TRUE)
            tmat2 = t(apply(tmat2, 1, cumsum))
            tmat2 = tmat2[indIds[ind] + Qk, ]
            Qk = max.col(runif(n) < tmat2, "first")
            speed = rnorm(n, object$mean[cbind(id[ind], Qk)], object$sd[cbind(id[ind], Qk)])
            tt = tt + len[k] * exp(-speed - E)
        }
    tt
}

