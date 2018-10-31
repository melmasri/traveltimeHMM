

traveltimeHMM <- function(logspeeds, trips, timeBins, linkIds, nQ = 1, model = c("HMM", 
    "trip-HMM", "no-dependence"), tolErr = 1, L = 10, maxItra = NULL, verbose = TRUE) {
    ## ################################################## Testing requirements
    if (length(logspeeds) != length(trips) || length(logspeeds) != length(timeBins)) 
        stop("Variabels logspeds, trips and timeBins are not equal in length!")
    
    if (!is.factor(trips)) 
        tripsFactor = factor(trips) else tripsFactor = trips
    
    if (!is.factor(linkIds)) 
        linkIds = factor(linkIds)
    
    ## ################################################## warnings for user
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

    nB = length(unique(timeBins))  # time bins
    nQ2 = nQ^2
    nObs = length(logspeeds)
    nlinks = nlevels(linkIds)  # number of links
    nTrips = nlevels(tripsFactor)
    
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
    initialObsIds = seq_along(tripsFactor)[!duplicated(tripsFactor)]  ## indices
    nInitalObsPerFactor = sapply(split(1:nlevels(tripsFactor), linkTimeFactor[initialObsIds]), 
        length)
    indexInitialObsLessMinObs = NULL
    
    ## finding the number of factors (link x timeBin) < L initial observations
    initialObsLessMinObs = which(nInitalObsPerFactor < L)
    indexInitialObsLessMinObs = sapply(gsub("[0-9.]+", "", names(initialObsLessMinObs)), 
        function(r) which(r == levels(timeFactor)))
    
    ## finding the number of factors (link x timeBin) that have only initial states,
    ## i.e cannot compute P(state_{k-1}, state_k | Obs})
    factorsWithOnlyInitialObs = which(sapply(split(1:nObs, linkTimeFactor), length) == 
        nInitalObsPerFactor)
    indexFactorsWithOnlyInitialObs = sapply(gsub("[0-9.]+", "", names(factorsWithOnlyInitialObs)), 
        function(r) which(r == levels(timeFactor)))
    
    ## ################################################## Storage variables Gaussian
    ## parameter matrices
    meanGauss = matrix(1:nQ - 1, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    varGauss = matrix(1:nQ, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    
    ## markov transition matrices and initial states per link
    tmat = matrix(rep(rep(1, nQ)/nQ, each = nQ), nrow = nB * nlinks, ncol = nQ2, 
        byrow = TRUE)
    init = matrix(rep(1, nQ)/nQ, nrow = nB * nlinks, ncol = nQ, byrow = TRUE)
    
    
    ## Empyt variables
    probStates <- tmatNew <- initNew <- meanGaussNew <- varGaussNew <- NULL
    tmatAverage <- initAverage <- meanSigAverage <- NULL
    
    error <- function(A, B) if (is.null(A) || is.null(B)) 
        0 else sqrt(sum(abs(A - B)^2))
    ## ################################################## Start of simulation
    ## Parameter estimation
    iter = 0  # iteration count
    tstart = Sys.time()  # starting time
    print(paste("model is :", model))
    repeat {
        ## forward-backward
        if (grepl("HMM", model)) {
            probTran = tmat[paramId, ] * pmax(dnorm(logspeeds,
                meanGauss[paramId,], varGauss[paramId, ]), 0.001)[, rep(1:nQ, nQ)]
            fb = tapply(1:nObs, tripsFactor, function(r)
                forwardback(probTran[r,], init[paramId[r[1]], ]))
            
            ## probability of each state per observation
            alpha = matrix(unlist(lapply(fb, function(r) r$alpha), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # forwad prob. normalized
            beta = matrix(unlist(lapply(fb, function(r) r$beta), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # backward prob. normalized
            probStates = normalizeR(alpha * beta + 1e-05, m = nObs, n = nQ)  # smoothed and normalized prob
            
            ## calculating initial state probability
            initNew = initial_est(probStates, linkTimeFactor, initialObsIds)
            if (!is.null(indexInitialObsLessMinObs)) {
                initAverage = initial_est(probStates, timeFactor, initialObsIds)
                initNew[initialObsLessMinObs, ] = initAverage[indexInitialObsLessMinObs, 
                  ]
            }
            ## calculating conditional joint transition probability P(state_k, state_{k-1} |
            ## Obs)
            probJointStates = normalizeR(alpha[-nObs, rep(1:nQ, each = nQ)] * (probTran * 
                beta[, rep(1:nQ, nQ)])[-1, ], m = nObs - 1, nQ2)
            probJointStates[initialObsIds[-1] - 1, ] <- 0  # setting the overlapping multiplication between trips to 0
            
            numGamma = lapply(split(probJointStates, linkTimeFactor[-1]), function(r) .colSums(r, 
                length(r)/nQ2, nQ2))
            numGamma = matrix(unlist(numGamma, use.names = FALSE), ncol = nQ2, byrow = TRUE)
            denGamma = probStates
            denGamma[initialObsIds[-1] - 1, ] <- 0
            denGamma = lapply(split(denGamma[-nObs, ], linkTimeFactor[-1]), function(r) .colSums(r, 
                length(r)/nQ, nQ))
            denGamma = matrix(unlist(denGamma, use.names = FALSE), ncol = nQ, byrow = TRUE)
            tmatNew = numGamma/denGamma[, rep(1:nQ, each = nQ)]
            
            if (!is.null(indexInitialObsLessMinObs)) {
                numGamma = lapply(split(probJointStates, timeFactor[-1]), function(r) .colSums(r, 
                  length(r)/nQ2, nQ2))
                numGamma = matrix(unlist(numGamma, use.names = FALSE), ncol = nQ2, 
                  byrow = TRUE)
                denGamma = probStates
                denGamma[initialObsIds[-1] - 1, ] <- 0
                denGamma = lapply(split(denGamma[-nObs, ], timeFactor[-1]), function(r) .colSums(r, 
                  length(r)/nQ, nQ))
                denGamma = matrix(unlist(denGamma, use.names = FALSE), ncol = nQ, 
                  byrow = TRUE)
                tmatAverage = numGamma/denGamma[, rep(1:nQ, each = nQ)]
                tmatNew[initialObsLessMinObs, ] <- tmatAverage[indexInitialObsLessMinObs, 
                  ]  # applying to all factors with less than min Obs
                if (length(factorsWithOnlyInitialObs)) {
                  tmatNew[factorsWithOnlyInitialObs, ] <- tmatAverage[indexFactorsWithOnlyInitialObs, 
                    ]
                }
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
        meanGaussNew = meanSig$mean
        varGaussNew = meanSig$sigma
        ord = order_states(meanSig$mean)
        if (!is.null(ord$order)) {
            meanGaussNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) meanGaussNew[ord$toSort[r], 
                ord$order[r, ]], USE.NAMES = FALSE))
            varGaussNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) varGaussNew[ord$toSort[r], 
                ord$order[r, ]], USE.NAMES = FALSE))
        }
        
        ## Calculating || ThetaNew - Theta ||
        maxError = sum(error(init, initNew), error(tmat, tmatNew), error(meanGauss, 
            meanGaussNew), error(varGaussNew, varGauss))
        ## re-positioning parameters
        tmat = tmatNew
        init = initNew
        meanGauss = meanGaussNew
        varGauss = varGaussNew
        iter = iter + 1
        print(paste(round(maxError, 3), "total error in iteration", iter, " @", format(Sys.time() - 
            tstart)))
        
        ## breaking loop on convergence
        if (maxError <= tolErr) {
            print(paste("Parameters converged at iteration ", iter))
            break
        }
        if (!is.null(maxItra) & iter >= maxItra) {
            print("Reached maximum number of iterations")
            break
        }
    }
    ## naming rows of variabels
    ## rownames(tmat)<-rownames(init)<-rownames(varGauss)<-rownames(meanGauss)<-
    ## levels(linkTimeFactor) rownames(tmatAverage)<-rownames(initAverage)<-rownames(
    ## meanSigAverage$mean)<-rownames(meanSigAverage$sigma)<-levels(timeFactor)
    
    ## returning variabels
    invisible(list(factors = linkTimeFactor, tmat = tmat, init = init, sd = sqrt(varGauss), 
        mean = meanGauss, nQ = nQ, nB = nB, nObs = nObs, model = model, tmatAverage = tmatAverage, 
        initAverage = initAverage, meanAverage = meanSigAverage$mean, sdAverage = sqrt(meanSigAverage$sigma)))
}
