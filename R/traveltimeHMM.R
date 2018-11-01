#' Commputes initial state probabilites 
#' \code{initial_est} return the state probabilities for each level in the passed \by{by_factor}.
#'
#' @param logspeeds 
#' @param trips
#' @param timeBins
#' @param linkIds
#' @param nQ Integer number of states
#' @param model 
#' @param max_iter
#' @param L
#' @param verbose
#'
#' @return 
#'
#' @examples
#' \dontrun{
#' }

traveltimeHMM <- function(logspeeds, trips, timeBins, linkIds, nQ = 1L, model = c("HMM", 
    "trip-HMM", "no-dependence"), tolErr = 1, L = 10, max_iter = NULL, verbose = TRUE) {

    ## #-------------------------------------------------- Testing requirements
    if (length(logspeeds) != length(trips) || length(logspeeds) != length(timeBins)) 
        stop("Variabels logspeds, trips and timeBins are not equal in length!")
    
    if (!is.factor(trips)) 
        trips = factor(trips) 
    
    if (!is.factor(linkIds)) 
        linkIds = factor(linkIds)
    
    ## #-------------------------------------------------- warnings for user
    maxSpeed = 130  # max speed km/h
    aux = which(logspeeds > log(maxSpeed/3.6))
    if (length(aux) > 0) 
        warning(paste("Many observations are higher than speed limit (130km/h)!, about", 
            paste0(round(100 * length(aux)/length(logspeeds), 2), "% of observations."), 
            " It is advised to remove these observations"))
    
    ## ################################################## Setup
    model <- match.arg(model)
    if (grepl("HMM", model) & nQ <= 1) 
        stop("Cannot use Hidden Markov Model with < 2 states!")

    nB <- length(unique(timeBins))  # time bins
    nQ2 <- nQ^2
    nObs <- length(logspeeds)
    nlinks <- nlevels(linkIds)  # number of links
    nTrips <- nlevels(trips)
    
    ## ################################################## Factors and sets
    timeFactor = interaction(timeBins, lex.order = TRUE)
    linkTimeFactor = interaction(linkIds, timeBins, lex.order = TRUE)  # factor of link id and time bin, ordered by lex.
    paramId = as.numeric(linkTimeFactor)
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
    index_only_init = sapply(gsub("[0-9.]+", "", names(only_init)),
        function(r) which(r == levels(timeFactor)))
    
    ## ################################################## Storage variables Gaussian
    ## parameter matrices
    mu_speed = matrix(1:nQ - 1, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    var_speed = matrix(1:nQ, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    
    ## markov transition matrices and initial states per link
    tmat = matrix(rep(rep(1, nQ)/nQ, each = nQ), nrow = nB * nlinks, ncol = nQ2, byrow = TRUE)
    init = matrix(rep(1, nQ)/nQ, nrow = nB * nlinks, ncol = nQ, byrow = TRUE)
    
    ## Empyt variables
    probStates <- mu_speedNew <- var_speedNew <- meanSigAverage <- NULL
    
    error <- function(A, B) if (is.null(A) || is.null(B)) 0 else sum(abs(A - B))

    ## #--------------------------------------------------
    ## Start of simulation

    ## Parameter estimation
    iter = 0  # iteration count
    tstart = Sys.time()  # starting time
    print(paste("model is :", model))
    repeat {
        ## forward-backward
        if (grepl("HMM", model)) {
            probTran = tmat[paramId, ] * pmax(dnorm(logspeeds,
                mu_speed[paramId,], var_speed[paramId, ]), 0.001)[, rep(1:nQ, nQ)]
            fb = tapply(1:nObs, trips, function(r)
                forwardback(probTran[r,], init[paramId[r[1]], ]))
            
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
            if (!is.null(index_init_L)) {
                tmat_impute = tmat_est(probJointStates, probStates, init_ids, timeFactor)
                tmatNew[init_L, ] <- tmat_impute[index_init_L, ]  # applying to all factors with less than min Obs
                if (length(only_init))
                    tmatNew[only_init, ] <- tmat_impute[index_only_init, ]
            }
        }
        
        ## calculating mean and variance of gaussian
        meanSig = gaussian_param_by_factor(logspeeds, linkTimeFactor, probStates)
        
        ## getting states with less than L factors
        if (!is.null(indexLinksLessMinObs)) {
            meanSigAverage = gaussian_param_by_factor(logspeeds, timeFactor, probStates)
            meanSig$mean[linksLessMinObs, ] = meanSigAverage$mean[indexLinksLessMinObs, 
                ]
            meanSig$sigma[linksLessMinObs, ] = meanSigAverage$sigma[indexLinksLessMinObs, 
                ]
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
        
        ## Calculating || ThetaNew - Theta ||
        iter_error = error(c(drop(init), drop(tmat), drop(mu_speed), drop(var_speed)),
            c(drop(initNew), drop(tmatNew), drop(mu_speedNew), drop(var_speedNew)))
        
        ## re-positioning parameters
        tmat <- tmatNew
        init <- initNew
        mu_speed <- mu_speedNew
        var_speed <- var_speedNew
        iter <- iter + 1
        
        print(paste(round(iter_error, 3),
                    "error in iteration", iter, " @",
                    format(Sys.time() - tstart)))
        
        ## breaking loop on convergence
        if (iter_error <= tolErr) {
            print(paste("Parameters converged at iteration ", iter-1))
            break
        }
        if (!is.null(max_iter) && iter >= max_iter) {
            print("Reached maximum number of iterations")
            break
        }
    }
    
    ## returning variabels
    invisible(list(factors = linkTimeFactor,
                   tmat = tmat,
                   init = init,
                   sd = sqrt(var_speed),
                   mean = mu_speed,
                   nQ = nQ,
                   nB = nB,
                   nObs = nObs,
                   model = model))
}
