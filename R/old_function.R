predict.traveltime.no_dependence <- function(object, linkIds, len, starttime, n = 1000) {
    if (object$model != "no-dependence") 
        stop("Object model is not no-dependence")
    fact = paste(linkIds[1], time_bins(starttime), sep = ".")
    id = which(levels(object$factors) == fact)
    speed = rnorm(n, object$mean[id, ], object$sd[id, ])
    tt = len[1] * exp(-speed)
    for (k in 2:length(linkIds)) {
        fact = paste(linkIds[k], time_bins(starttime + tt), sep = ".")
        if (length(unique(fact)) == 1) {
            id = which(levels(object$factors) == fact[1])
        } else {
            id = sapply(fact, function(s) which(levels(object$factors) == s))
        }
        speed = rnorm(n, object$mean[id, ], object$sd[id, ])
        tt = tt + len[k] * exp(-speed)
    }
    tt
}

predict.traveltime.HMM <- function(object, linkIds, len, starttime, n = 1000) {
    nQ = object$nQ
    fact = paste(linkIds[1], time_bins(starttime), sep = ".")
    id = which(levels(object$factors) == fact, useNames = FALSE)
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    tt = len[1] * exp(-speed)
    if (length(linkIds) > 1) 
        for (k in 2:length(linkIds)) {
            fact = paste(linkIds[k], time_bins(starttime + tt), sep = ".")
            if (length(unique(fact)) == 1) {
                id = which(levels(object$factors) == fact[1], useNames = FALSE)
                Qk = 1 * (runif(n) > matrix(object$tmat[id, ], ncol = nQ, byrow = TRUE)[Qk, 
                  1]) + 1
            } else {
                id = sapply(unique(fact), function(s) which(levels(object$factors) == 
                  s, useNames = FALSE), USE.NAMES = FALSE)
                id = id[as.numeric(as.factor(fact))]
                Qk = sapply(1:n, function(s) 1 * (runif(1) > matrix(object$tmat[id[s], 
                  ], ncol = nQ, byrow = TRUE)[Qk[s], 1]) + 1, USE.NAMES = FALSE)
            }
            speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
            tt = tt + len[k] * exp(-speed)
        }
    tt
}


predict.traveltime.HMM <- function(object, linkIds, len, starttime, n = 1000) {
    nQ = object$nQ
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
        sep = "."), useNames = FALSE)
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    tt = len[1] * exp(-speed)
    if (length(linkIds) > 1) 
        for (k in 2:length(linkIds)) {
            # this saves about 50ms for 117 routes (.4 ms per route)
            fact = paste(linkIds[k], time_bins(starttime + tt), sep = ".")
            id = sapply(unique(fact), function(s) which(levels(object$factors) == 
                s, useNames = FALSE), USE.NAMES = FALSE)
            ind = as.numeric(as.factor(fact))
            indIds = (1:length(id) - 1) * nQ
            ## creating tmat takes about 10ms for 2 rows or the full matrix to speed up create
            ## tmat2 before the loop or post-estimation for the whole Quebec dataset, for 1294
            ## obs it takes about an extra 11sec if at the top of the loop
            tmat2 = matrix(drop(t(object$tmat[id, ])), ncol = nQ, byrow = TRUE)
            tmat2 = t(apply(tmat2, 1, cumsum))
            tmat2 = tmat2[indIds[ind] + Qk, ]
            Qk = max.col(runif(n) < tmat2, "first")
            speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
            tt = tt + len[k] * exp(-speed)
        }
    tt
}

predict.traveltime.HMM.2 <- function(object, linkIds, len, starttime, n = 1000) {
    nQ = object$nQ
    lenIds = length(linkIds)
    ## uniIds = unlist(sapply(linkIds, function(r) 1:nB + nB*(r-1), simplify=FALSE,
    ## USE.NAMES=FALSE)) uniIds = rep(linkIds, 5) aux =
    ## strsplit(levels(object$factors), '[.]') uniIds =
    ## unlist(sapply(paste0('^',linkIds, '\\.'), function(r) grep(r,
    ## levels(object$factors)), USE.NAMES=FALSE)) uniFactors =
    ## levels(object$factors)[uniIds] levels(object$factors) aux =
    ## lapply(strsplit(levels(object$factors), '[.]'), function(r) r[1])
    ttQ_ <- function(k, travelTime, Q, linkIds, starttime) {
        if (k > lenIds) 
            return(travelTime)
        fact = paste(linkIds[k], time_bins(starttime + travelTime), sep = ".")
        ## id = match(fact, levels(object$factors))
        if (length(unique(fact)) == 1) {
            id = which(levels(object$factors) == fact[1], useNames = FALSE)
            ## id = uniIds[which(uniFactors == fact[1], useNames = FALSE)]
            Qin = 1 * (runif(n) > matrix(object$tmat[id, ], ncol = nQ, byrow = TRUE)[Q, 
                1]) + 1
        } else {
            id = sapply(unique(fact), function(s) which(levels(object$factors) == 
                s, useNames = FALSE), USE.NAMES = FALSE)
            id = id[as.numeric(as.factor(fact))]
            Qin = sapply(1:n, function(s) 1 * (runif(1) > matrix(object$tmat[id[s], 
                ], ncol = nQ, byrow = TRUE)[Q[s], 1]) + 1, USE.NAMES = FALSE)
        }
        speed = rnorm(n, object$mean[id, Qin], object$sd[id, Qin])
        ttQ_(k + 1, travelTime + len[k] * exp(-speed), Qin, linkIds, starttime)
    }
    
    ## first initial
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
        sep = "."), useNames = FALSE)
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    ttQ_(2, len[1] * exp(-speed), Qk, linkIds, starttime)
}


