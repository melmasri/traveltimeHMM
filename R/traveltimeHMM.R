#' \code{traveltimeHMM} estimates trips and link specific speed parameters from observed average speeds per unique trip and link.
#'
#' @param speeds A numeric vector of speed observations on the log-scale.
#' @param trips An integer or character vector of trip ids for each observation of \code{speed}.
#' @param timeBins A character vector of time bins for each observation of \code{speed}.
#' @param linkIds A vector of link ids (route or way) for each observation of \code{speed}.
#' @param nQ Integer number of states, default is \code{1}.
#' @param model Type of model as string, \code{trip-HMM} to use a hidden Markov model (HMM) with trip effect, \code{HMM} (default) is an HMM without trip effect,
#' \code{trip} is trip effect model without HMM, and \code{no-dependence} is model with link specific parameter only without an HMM nor a trip effect.
#' @param max.it An Integer for the maximum number of iterations to run for, default = \code{NULL}.
#' @param L An integer minimum number of observations per factor (\code{linkIds x timeBins}) to estimate the paramter for, \code{default = 10}. Factors that have less observations than \code{L} their estimates are imputed by the average over timeBins.
#' @param tol.err A numeric variable representing the level of tolerable distance between parameter estimates from consecutive iterations.
#' @param ... extra specific parameters, see details
#' @details 
#'
#'
#' 
#' @return  NULL
#'
#' @examples
#' \dontrun{
#' }
traveltimeHMM <- function(speeds, trips, timeBins, linkIds, nQ = 1L, model = c("HMM", 
    "trip-HMM","trip","no-dependence"), tol.err = 10, L = 10L, max.it = NULL, ...) {

    ## #-------------------------------------------------- Testing requirements
    if (length(speeds) != length(trips) || length(speeds) != length(timeBins)) 
        stop("Variabels logspeds, trips and timeBins are not equal in length!")
    
    ## #-------------------------------------------------- warnings for user
    maxSpeed = 130  # max speed km/h
    aux = which(speeds > log(maxSpeed/3.6))
    if (length(aux) > 0) 
        warning(paste("Many observations are higher than speed limit (130km/h)!, about", 
            paste0(round(100 * length(aux)/length(speeds), 2), "% of observations."), 
            " It is advised to remove these observations"))
    
    model <- match.arg(model)
    param <- list(...)
    if (grepl("HMM", model) & nQ <= 1) 
        stop("Cannot use Hidden Markov Model with < 2 states!")
    if(!grepl('HMM', model) & nQ > 1){
        warning('Cannot use nQ > 1 without an HMM, resetting nQ = 1', immediate.=TRUE, call.=FALSE)
        nQ =1
    }

    ## #-------------------------------------------------- setup
    if (!is.factor(trips)) 
        trips = factor(trips) 
    
    if (!is.factor(linkIds)) 
        linkIds = factor(linkIds)
    
    nB <- length(unique(timeBins))  # time bins
    nQ2 <- nQ^2
    nObs <- length(speeds)
    nlinks <- nlevels(linkIds)  # number of links
    nTrips <- nlevels(trips)
    
    ## ################################################## Factors and sets
    timeFactor = interaction(timeBins, lex.order = TRUE)
    linkTimeFactor = interaction(linkIds, timeBins, lex.order = TRUE)  # factor of link id and time bin, ordered by lex.
    obsId = as.numeric(linkTimeFactor)
    tripId = as.numeric(trips)
    nFactors = sapply(split(1:nObs, linkTimeFactor), length)
    
    ## which links (factors) with less than L
    indexLinksLessMinObs = NULL
    linksLessMinObs = which(nFactors < L)
    indexLinksLessMinObs = sapply(gsub("[0-9.]+", "", names(linksLessMinObs)), function(r) which(r == 
        levels(timeFactor)))
    
    ## finding the number of initial state observations per link x timeBin
    init_ids <- seq_along(trips)[!duplicated(trips)]  ## indices
    count_init <- sapply(split(1:nlevels(trips), linkTimeFactor[init_ids]),length)
    
    ## finding the number of factors (link x timeBin) < L initial observations
    init_L <- which(count_init < L)      
    index_init_L <- sapply(gsub("[0-9.]+", "", names(init_L)),function(r) which(r == levels(timeFactor)))
    
    ## finding the number of factors (link x timeBin) that have only initial states,
    ## i.e cannot compute P(state_{k-1}, state_k | Obs})
    only_init = which(sapply(split(1:nObs, linkTimeFactor), length) == count_init)
    index_only_init = sapply(gsub("[0-9.]+", "", names(only_init)),function(r) which(r == levels(timeFactor)))
    
    ## #-------------------------------------------------- variable holders
    ## speed variables
    mu_speed = matrix(1:nQ - 1, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    var_speed = matrix(1:nQ, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    
    ## Markov transition matrices and initial states per link
    ## default initial values
    init.0 <- rep(1,nQ)/nQ
    tmat.0 <-rep(rep(1, nQ)/nQ, nQ)
    ## sampling initial values if seed is passed
    if(!is.null(param$seed) && is.numeric(param$seed)){
        init.0 <- runif(nQ)
        init.0 <- init.0/sum(init.0)
        tmat.0 <- matrix(runif(nQ^2), nQ, nQ)
        tmat.0 <- c(t(tmat.0/rowSums(tmat.0)))
    }
    ## setting initial values if optional parameters are passed
    if(!is.null(param$tmat)){
        if(is.matrix(param$tmat) && dim(param$tmat) == c(nQ, nQ) && rowSums(param$tmat) == c(1,1)){
            tmat.0 <- c(t(param$tmat))
        }else
            warning('Initial transition values are not used, tmat must be nQ x nQ with rows summing to 1', immediate.=TRUE, call.=FALSE)
    }
    if(!is.null(param$init)){
        if(length(param$init) == nQ && sum(param$init) == 1){
            init.0 <- param$init
        }else
            warning('Initial state values are not used, init must be of length nQ and sum to 1', immediate.=TRUE, call.=FALSE)
    }
    ## creating initial matrices from inital values
    init <- matrix(init.0, nrow = nB * nlinks, ncol = nQ, byrow = TRUE)
    tmat <- matrix(tmat.0, nrow = nB * nlinks, ncol = nQ2, byrow = TRUE)
    
    ## Trip-effect 
    tau2 = 1
    if(grepl('trip', model)) E <- rnorm(nTrips, 0, sqrt(tau2)) else E <- numeric(nTrips)
    ## Empyt variables
    initNew <- tmatNew <- mu_speedNew <- var_speedNew <- NULL
    E_new <- E
    probStates <- NULL
    ## Error function
    error <- function(A, B) if (is.null(A) || is.null(B)) 0 else sum(abs(A - B))
    ## #--------------------------------------------------
    ## Start of simulation

    ## Parameter estimation
    iter = 0  # iteration count
    tstart = Sys.time()  # starting time
    cat("Running model", model, fill=TRUE)
    repeat {
        ## forward-backward
        if (grepl("HMM", model)) {
            probTran = tmat[obsId, ] * pmax(dnorm(speeds, mu_speed[obsId,], var_speed[obsId, ]), 0.001)[, rep(1:nQ, nQ)]
            fb = tapply(1:nObs, trips, function(r)
                forwardback(probTran[r,], init[obsId[r[1]], ]))
            
            ## probability of each state per observation
            alpha = matrix(unlist(lapply(fb, function(r) r$alpha), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # forwad prob. normalized
            beta = matrix(unlist(lapply(fb, function(r) r$beta), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # backward prob. normalized
            probStates = normalizeR(alpha * beta + 1e-05, m = nObs, n = nQ)  # smoothed and normalized prob
            
            ## calculating initial state probability
            initNew = initial_est(probStates, linkTimeFactor, init_ids)
            if (length(index_init_L)) {
                init_impute = initial_est(probStates, timeFactor, init_ids)
                initNew[init_L, ] = init_impute[index_init_L, ]
            }

            ## calculating conditional joint transition probability P(state_k, state_{k-1} |
            ## Obs)
            probJointStates = normalizeR(alpha[-nObs, rep(1:nQ, each = nQ)] *
                                             (probTran * beta[, rep(1:nQ, nQ)])[-1, ], m = nObs - 1, nQ2)
            probJointStates[init_ids[-1] - 1, ] <- 0  # setting the overlapping multiplication between trips to 0

            ## computing transition matrix estimate and imputing < threshold with second factor
            tmatNew = tmat_est(probJointStates, probStates, init_ids, linkTimeFactor)
            ## if (!is.null(index_init_L)) {
            ##     tmat_impute = tmat_est(probJointStates, probStates, init_ids, timeFactor)
            ##     tmatNew[init_L, ] <- tmat_impute[index_init_L, ]  # applying to all factors with less than min Obs
            ## }
            if (!is.null(indexLinksLessMinObs)) {
                tmat_impute = tmat_est(probJointStates, probStates, init_ids, timeFactor)
                tmatNew[linksLessMinObs, ] <- tmat_impute[indexLinksLessMinObs ]  # applying to all factors with less than min Obs
                if (length(only_init))
                    tmatNew[only_init, ] <- tmat_impute[index_only_init, ]
            }
            
        }
        
        ## calculating mean and variance of gaussian
        meanSig = gaussian_param_by_factor(speeds - E[tripId], linkTimeFactor, probStates)
        ## getting states with less than L factors
        if (!is.null(indexLinksLessMinObs)) {
            meanSigAverage = gaussian_param_by_factor(speeds - E[tripId], timeFactor, probStates)
            meanSig$mean[linksLessMinObs, ] = meanSigAverage$mean[indexLinksLessMinObs, ]
            meanSig$sigma[linksLessMinObs, ] = meanSigAverage$sigma[indexLinksLessMinObs, ]
        }
        
        ## sorting mean of Gaussian parameter to make the HMM states identifiable
        mu_speedNew = meanSig$mean
        var_speedNew = meanSig$sigma
        ord = order_states(meanSig$mean)
        if (!is.null(ord$order)) {
            mu_speedNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) mu_speedNew[ord$toSort[r], 
                ord$order[r, ]], USE.NAMES = FALSE))
            var_speedNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) var_speedNew[ord$toSort[r], 
                ord$order[r, ]], USE.NAMES = FALSE))
        }
        

        ## Calculating E-effect (trip specific effect parameters)
        if(grepl('trip', model)){
            ## Calculating E-trip variance
            tau2 <- .colSums(E^2, m = nTrips, n=1)/nTrips

            ## Updaing E-trip per trip
            dummy <- if(is.null(probStates)) 1/var_speedNew[obsId,] else  probStates/var_speedNew[obsId,]
            a <- .rowSums(dummy, m = nObs, n = nQ)
            h <- .rowSums(dummy * mu_speedNew[obsId,], m = nObs, n = nQ)
            dummy1 <- vapply(split(speeds * a - h, trips), function(r) .colSums(r, m = length(r), n=1), numeric(1), USE.NAMES = FALSE)
            dummy2 <- vapply(split(a, trips), function(r) .colSums(r, m = length(r), n=1), numeric(1), USE.NAMES = FALSE)
            E_new <-  dummy1 / (1/tau2 + dummy2)
        }
        
        ## Calculating || ThetaNew - Theta ||
        A = c(mu_speed[-linksLessMinObs,], var_speed[-linksLessMinObs,])
        B = c(mu_speedNew[-linksLessMinObs,], var_speedNew[-linksLessMinObs,])
        if(grepl('HMM', model)){
            A = c(A, init[-init_L,], tmat[-unique(linksLessMinObs,only_init),])
            B = c(B, initNew[-init_L,], tmatNew[-unique(linksLessMinObs,only_init),])
        }
        if(grepl('trip',model)) { A = c(A, E); B = c(B,E_new)}
        iter_error = error(A, B)

        ## re-positioning parameters
        tmat <- tmatNew
        init <- initNew
        mu_speed <- mu_speedNew
        var_speed <- var_speedNew
        E <- E_new
        iter <- iter + 1

        if(!is.null(max.it) && iter==1)
            cat('Expected completion of', max.it, 'iterations in', format((max.it-1) * (Sys.time() - tstart), digits = 3), fill=TRUE)

        cat(round(iter_error,2), "error in iteration", iter, "@", format(Sys.time() - tstart, digits = 3), fill=TRUE)

                
        ## breaking loop on convergence
        if (iter_error <= tol.err) {
            cat("Parameters converged at iteration" , iter-1, fill=TRUE)
            break
        }
        if (!is.null(max.it) && iter >= max.it) {
            cat("Reached maximum number of iterations", fill=TRUE)
            break
        }
    }
    
    ## returning variabels
    invisible(list(factors = linkTimeFactor,
                   tmat = tmat,
                   init = init,
                   sd = sqrt(var_speed),
                   mean = mu_speed,
                   tau  = sqrt(tau2),
                   E = E,
                   nQ = nQ,
                   nB = nB,
                   nObs = nObs,
                   model = model))
}
