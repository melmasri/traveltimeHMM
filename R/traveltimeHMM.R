#' \code{traveltimeHMM} estimates trips and link specific speed parameters from observed average speeds per unique trip and link.
#'
#' @param logspeeds A numeric vector of speed observations (in km/h) on the (natural) log-scale.
#' @param trips An integer or character vector of trip ids for each observation of \code{speed}.
#' @param timeBins A character vector of time bins for each observation of \code{speed}.
#' @param linkIds A vector of link ids (route or way) for each observation of \code{speed}.
#' @param nQ Integer number of states, default is \code{1}.
#' @param model Type of model as string, \code{trip-HMM} to use a hidden Markov model (HMM) with trip effect, \code{HMM} (default) is an HMM without trip effect,
#' \code{trip} is trip effect model without HMM, and \code{no-dependence} is model with link specific parameter only without an HMM nor a trip effect.
#' @param tol.err A numeric variable representing the level of tolerable distance between parameter estimates from consecutive iterations.
#' @param L An integer minimum number of observations per factor (\code{linkIds x timeBins}) to estimate the parameter for, \code{default = 10}. Factors that have less observations than \code{L} their estimates are imputed by the average over timeBins.
#' @param max.it An integer for the maximum number of iterations to run for, default = \code{NULL}.
#' @param maxSpeed A float for the maximum speed in km/h, on the linear scale (not the log-scale, unlike for \code{logspeeds}).
#' @param ... extra specific parameters, see details
#' @details NULL
#' 
#' @return \code{traveltimeHMM} returns a list of the following parameters
#' 
#' \item{factors}{a factor of interactions (linkId x timeBin) with the same length of observations, and with levels corresponding to unique factors}
#' \item{trip}{a factor of trips.}
#' \item{tmat}{a transition matrix with rows corresponding to \code{levels(factors)}, with columns being the row-wise transition matrix of that factor. For example, \code{matrix(tmat[1,], ncol = nQ, nrow = nQ, byrow = TRUE)} is the transition matrix of \code{levels(factors)[1]}.}
#' \item{init}{a initial state probability matrix with rows corresponding to \code{levels(factors)}, and columns to the \code{nQ} states.}
#' \item{sd}{a matrix of standard deviations estimates with rows corresponding to \code{levels(factors)}, and columns to standard deviation estimates of the \code{nQ} states.}
#' \item{mean}{a matrix of mean estimates with rows corresponding to \code{levels(factors)}, and columns to mean estimates for the \code{nQ} states.}
#' \item{tau}{a numeric variable for the standard deviation estimate of the trip effect parameter \code{E}.}
#' \item{E}{a numeric vector of trip effect estimates corresponding to \code{levels(trip)}.}
#' \item{nQ}{an integer number of states.}
#' \item{nB}{an integer number of unique time bins.}
#' \item{nObs}{an integer number of observations.}
#' \item{model}{the type of model used.}
#'
#' @examples
#' \dontrun{
#' data(tripset)
#' ?traveltimeHMM  # for help
#' fit <- traveltimeHMM(tripset$logspeed, tripset$trip, tripset$timeBin, tripset$linkId, nQ = 2, max.it = 2)
#' single_trip <- subset(tripset, trip==2700)
#' pred <- predict.traveltime(fit, single_trip,single_trip$time[1])
#' hist(pred)      # histogram of prediction samples
#' mean(pred)      # travel time point estimate
#' sum(single_trip$traveltime)    # observed traveltime
#' }
#' @export
traveltimeHMM <- function(logspeeds, trips, timeBins, linkIds, nQ = 1L,
                          model = c("HMM", "trip-HMM","trip","no-dependence"),
                          tol.err = 10, L = 10L, max.it = 20,verbose = FALSE,
                          maxSpeed = NULL, seed = NULL, transition_matrix = NULL,
                          initial_state_values = NULL) {

    # SECTION A - Parameter validation and related processing
    
    # param <- list(...)  # Put this line back into code later only if required.
    # Ideally all parameters should be explicit in the interface.  EG 2019/05/29
  
    # A.1 Stop if 'logspeeds', 'trips' and 'timeBins' have different sizes
    if (length(logspeeds) != length(trips) || length(logspeeds) != length(timeBins)) 
        stop("Parameter vectors for 'logspeeds', 'trips' and 'timeBins' are not equal in length!")
  
    # A.2 Parameter 'model'
    # Get model from "HMM", "trip-HMM","trip" or "no-dependence"; stop if invalid
    model <- tryCatch(match.arg(model),error=function(cond){
      stop("Parameter 'model' should be one of “HMM”, “trip-HMM”, “trip”, “no-dependence”")
    })

    # Output chosen model
    if (verbose) message(paste('Using model ', model))
    
    # A.3 Parameter 'nQ'
    # Hidden Markov Model needs at least 2 states whilst other models require exactly 1 state.
    # All other values are forbidden.  In case of invalid value, we set nQ to 1 if permitted;
    # otherwise we stop.
    
    if (grepl("HMM", model) & nQ <= 1) 
        stop("Cannot use Hidden Markov Model with < 2 states!")
    if (nQ < 0) {
      warning("Cannot use nQ <1, resetting nQ = 1", immediate.=TRUE, call.=FALSE)
      nQ <- 1
    }
    if(!grepl('HMM', model) & nQ > 1){
        warning('Cannot use nQ > 1 without an HMM, resetting nQ = 1', immediate.=TRUE, call.=FALSE)
        nQ <- 1
    }
    
    # A.4 Parameter 'maxSpeed'
    # Set default value if not specified in call
    if(is.null(maxSpeed)){
        message('maxSpeed is not specified, setting at default value: 130 km/h')
        maxSpeed = 130
    }
    
    # We count occurrences of speeds in excess of maxSpeed.  If any, we then
    # report the percentage of such occurrences in a warning to the user.
    speedings = which(logspeeds > log(maxSpeed/3.6))
    if (length(speedings) > 0) 
        warning(paste("Many observations are higher than speed limit (130km/h)!, about", 
            paste0(round(100 * length(speedings)/length(logspeeds), 2), "% of observations."), 
            " It is advised to remove these observations."),  immediate.=TRUE, call.=FALSE)

    # SECTION B - Data preparation
    
    # Convert 'trips' and 'linkIds' to factors if required
    if (!is.factor(trips)) 
        trips = factor(trips) 
    if (!is.factor(linkIds)) 
        linkIds = factor(linkIds)
    
    nQ2 <- nQ^2 # Compute nQ2, the square of nQ, which will be used often
    nObs <- length(logspeeds) # nObs now contains the number of observations
    nlinks <- nlevels(linkIds) # nlinks now contains the number of distinct links
    nTrips <- nlevels(trips) # nTrips now contains the number of distinct trips
    
    # timeFactor = interaction(timeBins, lex.order = TRUE) # Replaced by next line
    timeFactor = factor(timeBins)
    linkTimeFactor = interaction(linkIds, timeBins, lex.order = TRUE)  # factor of link id and time bin, using lexicographic order.
    obsId = as.numeric(linkTimeFactor) # CAUTION with that approach.  ÉG 2019/05/30
    tripId = as.numeric(trips) # CAUTION with that approach.  ÉG 2019/05/30
    nFactors = sapply(split(1:nObs, linkTimeFactor), length)       # Number of observations per linkTimeFactor
    
    # Preparation of imputation variables - refers to eq. (6) in Woodard et al.
    # indexLinksLessMinObs = NULL # Useless line it seems, confirm removal.  ÉG 2019/06/03
    linksLessMinObs = which(nFactors < L) # vector of indices of all factors for which we need to impute 
    indexLinksLessMinObs = sapply(gsub("[0-9.]+", "", names(linksLessMinObs)), function(r) which(r == 
        levels(timeFactor))) # vector of all timeBins with counts of occurrences of the corresponding linkTimeFactor...
                             # We need to discuss this one.  ÉG 2019/06/03
    
    # Taking into account the initial state observation (initial occurrence) for each link for each time bin...
    init_ids <- seq_along(trips)[!duplicated(trips)]  # vector of indices of the initial state observation for each trip
    count_init <- sapply(split(1:nlevels(trips), linkTimeFactor[init_ids]),length) # vector of counts of initial
                                                                                   # state observations for each linkTimeFactor
    
    # finding the number of factors (link x timeBin) < L initial observations
    init_L <- which(count_init < L) # vector of indices of all linkTimeFactors with insufficient occurrences
    index_init_L <- sapply(gsub("[0-9.]+", "", names(init_L)),function(r) which(r == levels(timeFactor))) # vector of counts
                                                                                   # of occurrences for each timeBin
                                                                                   # below L occurrences.  To be discusse.
                                                                                   # ÉG 2019/06/03
    
    ## finding the number of factors (link x timeBin) that have only initial states,
    ## i.e cannot compute P(state_{k-1}, state_k | Obs})
    only_init = which(sapply(split(1:nObs, linkTimeFactor), length) == count_init) # get indices of linkTimeFactors for
                                                                                   # unique trips, i.e. those for which
                                                                                   # the number ofoccurrences equals
                                                                                   # the number of initial occurrences
    index_only_init = sapply(gsub("[0-9.]+", "", names(only_init)),function(r) which(r == levels(timeFactor))) # vector of
                                                                                   # counts of occurrences for each timeBin
                                                                                   # involved in unique trips only.
    
    # Travel speed variables:
    # mu_speed and var_speed are the mean and variance matrices in eq. (3) in Woodard et al.
    # Each matrix has nQ columns and as many rows as there are linkTimeFactors.
    # All rows are identical within each matrix:
    # # for mu_speed, column values are 0..nQ-1
    # # for var_speed, column values are 1..nQ
    mu_speed = matrix(1:nQ - 1, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    var_speed = matrix(1:nQ, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    
    # Markov transition matrices and initial states per link:
    # refer to eq. (4) in Woodard et al.
    # init refers to small 'gamma' whilst tmat refers to big 'gamma'.
    
    # default initial values init.0 and tmat.0
    # WARNING: those values might be superseded, see later.  To be discussed.  ÉG 2019/06/03
    init.0 <- rep(1,nQ)/nQ
    tmat.0 <-rep(rep(1, nQ)/nQ, nQ)
    
    # Here we supersede the values just set for init.0 and tmat.0 if seed is passed
    if(!is.null(seed) && is.numeric(seed)){
        init.0 <- runif(nQ) # gives a vector of random values between 0 and 1
        init.0 <- init.0/sum(init.0) # normalizes the vector
        tmat.0 <- matrix(runif(nQ^2), nQ, nQ) # gives a matrix of random values between 0 and 1
        tmat.0 <- c(t(tmat.0/rowSums(tmat.0))) # normalizes each row and transforms into vector
    }
    # Here we attempt to supersede the value for the transition matrix if one is passed
    if(!is.null(transition_matrix)){
        if(is.matrix(transition_matrix) && dim(transition_matrix) == c(nQ, nQ) && rowSums(transition_matrix) == rep(1,nQ)){
            tmat.0 <- c(t(transition_matrix))
        }else
            warning('Initial transition values are not used, tmat must be nQ x nQ with rows summing to 1', immediate.=TRUE, call.=FALSE)
    }
    
    # Here we attempt to supersede the value for the initial state vector if one is passed
    if(!is.null(initial_state_values)){
        if(length(initial_state_values) == nQ && sum(initial_state_values) == 1){
            init.0 <- initial_state_values
        }else
            warning('Initial state values are not used, init must be of length nQ and sum to 1', immediate.=TRUE, call.=FALSE)
    }
    # Creating initial matrices 'init' and 'tmat' from initial values 'init.0' and 'tmat.0'.
    # init has nQ columns whilst tmat has nQ*nQ columns.
    # Each matrix has as many rows as there are linkTimeFactors.
    # All rows are identical within each matrix.  Each row of init has the values of init.0
    # # whereas each row of tmat has the values of tmat.0
    
    init <- matrix(init.0, nrow = nB * nlinks, ncol = nQ, byrow = TRUE)
    tmat <- matrix(tmat.0, nrow = nB * nlinks, ncol = nQ2, byrow = TRUE)
    
    # Modelling the trip effect E: see eq. (1) from Woodard et al.
    tau2 = 1 # set initial value for tau2: see eq. (5) from Woodard et al.
    
    # We create the vector E of length nTrips as follows:
    # if model is "trip-HMM" or "trip" then each Ei has a random value of mean 0 and variance tau2;
    # otherwise each Ei has a value of zero.
    if(grepl('trip', model)) E <- rnorm(nTrips, 0, sqrt(tau2)) else E <- numeric(nTrips)
    
    ## Initialize variables - TO DOCUMENT FURTHER - ÉG 2019/06/03
    initNew <- tmatNew <- mu_speedNew <- var_speedNew <- probStates <- NULL
    E_new <- E
    iter = 0  # iteration count
    tstart = Sys.time()  # starting time
    
    ## Definition of error function as the absolute value of the difference between both parameters
    error.fun <- function(A, B) if (is.null(A) || is.null(B)) 0 else sum(abs(A - B))
    
    ###
    # Start of simulation
    # This corresponda to Algorithm 1 (which we call Algo1 hereafter) in Woodard et al.
    ###
    
    # Message to user:
    message('Model ', model, ' with ', nTrips, ' trips over ',
            nlinks, ' roads and ', nB, ' time bins...')
    
    repeat {
        if (grepl("HMM", model)) { 
            probTran = tmat[obsId, ] * pmax(dnorm(logspeeds, E[tripId] + mu_speed[obsId,], var_speed[obsId, ]), 0.001)[, rep(1:nQ, nQ)]
            fb = tapply(1:nObs, trips, function(r)
                forwardback(probTran[r,], init[obsId[r[1]], ]))
            
            ## probability of each state per observation
            alpha = matrix(unlist(lapply(fb, function(r) r$alpha), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # forward prob. normalized
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
            if (!is.null(indexLinksLessMinObs)) {
                tmat_impute = tmat_est(probJointStates, probStates, init_ids, timeFactor)
                tmatNew[linksLessMinObs, ] <- tmat_impute[indexLinksLessMinObs ]  # applying to all factors with less than min Obs
                if (length(only_init))
                    tmatNew[only_init, ] <- tmat_impute[index_only_init, ]
            }
            
        }
        
        ## calculating mean and variance of Gaussian
        meanSig = gaussian_param_by_factor(logspeeds - E[tripId], linkTimeFactor, probStates)
        ## getting states with less than L factors
        if (!is.null(indexLinksLessMinObs)) {
            meanSigAverage = gaussian_param_by_factor(logspeeds - E[tripId], timeFactor, probStates)
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

            ## Updating E-trip per trip
            dummy <- if(is.null(probStates)) 1/var_speedNew[obsId,] else  probStates/var_speedNew[obsId,]
            a <- .rowSums(dummy, m = nObs, n = nQ)
            h <- .rowSums(dummy * mu_speedNew[obsId,], m = nObs, n = nQ)
            dummy1 <- vapply(split(logspeeds * a - h, trips), function(r) .colSums(r, m = length(r), n=1), numeric(1), USE.NAMES = FALSE)
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
        iter_error = error.fun(A, B)

        ## re-positioning parameters
        tmat <- tmatNew
        init <- initNew
        mu_speed <- mu_speedNew
        var_speed <- var_speedNew
        E <- E_new
        iter <- iter + 1

        if(iter==1 && (!is.null(max.it) || verbose)){
            if(!is.null(max.it)){
                message('Expected completion of ', max.it, ' iterations in ',
                        format((max.it-1) * (Sys.time() - tstart), digits = 3))
            }else{
                message('Expected time per iteration ',
                        format((Sys.time() - tstart), digits = 3))
            }
        }
        if(verbose)
            message(round(iter_error,2), " error in iteration ", iter, " @ ", format(Sys.time() - tstart, digits = 3))

                
        ## breaking loop on convergence
        if (iter_error <= tol.err) {
            message("Parameters converged at iteration " , iter-1)
            break
        }
        if (!is.null(max.it) && iter >= max.it) {
            message("Reached maximum number of iterations")
            break
        }
    }
    
    ## returning variables
    obj <- list(factors = linkTimeFactor,
                trip = trips,
                tmat = tmat,
                init = init,
                sd = sqrt(var_speed),
                mean = mu_speed,
                tau  = sqrt(tau2),
                E = E,
                nQ = nQ,
                nB = nB,
                nObs = nObs,
                model = model)
    class(obj) <- append(class(obj),"traveltime")
    invisible(obj)
}

