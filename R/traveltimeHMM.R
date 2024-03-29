#' Estimate trip- and link- specific speed parameters from observed average speeds
#' 
#' \code{traveltimeHMM} estimates trip- and link- specific speed parameters from observed average speeds per unique trip and link.
#'
#' @param logspeeds A numeric vector of speed observations (in km/h) on the (natural) log-scale.  Needs
#'   to be provided if \code{data} is \code{NULL}, otherwise it will be overwritten by \code{logspeeds} column in \code{data}.  Default is \code{NULL}.
#' @param trips An integer or character vector of trip ids for each observation of \code{speed}.  Needs
#'   to be provided if \code{data} is \code{NULL}, otherwise it will be overwritten by \code{trips} column in \code{data}. Default is \code{NULL}.
#' @param timeBins A character vector of time bins for each observation of \code{speed}. Needs
#'   to be provided if \code{data} is \code{NULL}, otherwise it will be overwritten by \code{timeBins} column in \code{data}. Default is \code{NULL}.
#' @param linkIds A vector of link ids (route or way) for each observation of \code{speed}.  Needs
#'   to be provided if \code{data} is \code{NULL}, otherwise it will be overwritten by \code{linkIds} column in \code{data}. Default is \code{NULL}.
#' @param data A data frame or equivalent object that contains one column for each of
#'   \code{logspeeds}, \code{trips}, \code{timeBins} and \code{linkIds}.  Default is \code{NULL}.  This
#'   argument is mutually exclusive with the full joint set of \code{logspeeds}, \code{trips},
#'   \code{timeBins} and \code{linkIds}. One can either pass the latter vectors explicitly or as a dataframe \code{data}.
#' @param nQ An integer corresponding to the number of different congestion states that the traversal
#'   of a given link can take corresponding to \code{|{1...., Q}|}.  Hidden Markov models (see table 2)
#'   require \code{nQ >= 2} whilst other models require exactly \code{nQ = 1}.  As an example,
#'   Woodard uses \code{nQ = 1}.  Default is \code{1}.
#' @param model Type of model as string, \code{trip-HMM} to use a hidden Markov model (HMM) with trip effect, \code{HMM} (default) is an HMM without trip effect,
#'   \code{trip} is trip effect model without HMM, and \code{no-dependence} is model with link specific
#'   parameter only without an HMM nor a trip effect.
#' @param tol.err A numeric variable representing the threshold under which the estimation algorithm will
#'   consider it has reached acceptable estimate values.  Default is \code{10}.
#' @param L An integer minimum number of observations per factor (\code{linkIds x timeBins}) to estimate
#'   the parameter for.  Default is \code{10}. Factors that have fewer total observations or initial state
#'   observations than \code{L} have their estimates imputed using time bin data, without regard to road
#'   link data.
#' @param max.it An integer for the upper limit of the iterations to be performed by the estimation
#'   algorithm.  Default is \code{20}.
#' @param verbose A boolean that triggers verbose output.  Default is \code{FALSE}.
#' @param max.speed An optional float for the maximum speed in km/h, on the linear scale
#'   (not the log-scale, unlike for \code{logspeeds}).  Default is \code{NULL} which results in a maximum speed
#'   of 130 km/h.  If there is any occurrence of speed in excel of \code{max.speed}, then the system issues to the
#'   user a warning that includes the percentage of such occurrences.
#' @param seed An optional float for the seed used for the random generation of the first Markov transition matrix
#'   and initial state vector.  Default is \code{NULL}.  If not provided, then those objects are generated deterministically.
#'   The effect of \code{seed} is cancelled by tmat.p or init.p when provided.
#' @param tmat.p An optional starting value for the Markov transition matrix \eqn{\Gamma} of size \code{nQ x nQ}
#'   with rows summing to \code{1}.  Default is \code{NULL}
#' @param init.p An optional starting value for the Markov initial state vector \eqn{\gamma}
#'   of size \code{nQ} with elements summing to \code{1}.  Default is \code{NULL}.
#' @details NULL
#' 
#' @return \code{traveltimeHMM} returns a list of the following parameters.
#' \item{factors}{a factor of interactions (linkId x timeBin) of length \code{nObs}.  Factors are in the format \code{linkId.timeBin}.}
#' \item{trip}{a factor of trip IDs.}
#' \item{tmat}{a transition matrix with rows corresponding to \code{levels(factors)}, and with columns being the row-wise transition matrix of that factor. For example, \code{matrix(tmat[1,], ncol = nQ, nrow = nQ, byrow = TRUE)} is the transition matrix of \code{levels(factors)[1]}.  \code{NULL} if hidden Markov modelling is not handled by the selected model type.}
#' \item{init}{a initial state probability matrix with rows corresponding to \code{levels(factors)}, and columns to the \code{nQ} states.  \code{NULL} if hidden Markov modelling is not handled by the selected model type.}
#' \item{sd}{a matrix of standard deviations estimates for the natural logarithm of the speed (in km/h), with rows corresponding to \code{levels(factors)}, and columns to standard deviation estimates for the \code{nQ} states.}
#' \item{mean}{a matrix of mean estimates for the natural logarithm of the speed (in km/h), with rows corresponding to \code{levels(factors)}, and columns to mean estimates for the \code{nQ} states.}
#' \item{tau}{a numeric variable for the standard deviation estimate for the trip effect parameter \code{logE}.   Equals \code{1} if trip effect is not handled by the selected model type.}
#' \item{logE}{a numeric vector of trip effect estimates corresponding to \code{levels(trip)}.  Values are set to \code{0} if trip effect is not handled by the selected model type. Units are the same as for \code{logspeeds}.}
#' \item{nQ}{an integer corresponding to the number of different congestion states, equal to the parameter nQ that was passed in the function call.}
#' \item{nB}{an integer number of unique time bins.}
#' \item{nObs}{an integer number of observations.}
#' \item{model}{a character string corresponding to the type of model used.}
#'
#' @examples
#' \dontrun{
#' data(tripset)
#' 
#' # Fit an HMM model with 2 hidden congestion states and 20 algorithm iterations
#' fit <- traveltimeHMM(tripset$logspeed, tripset$tripID, tripset$timeBin, 
#'                       tripset$linkID, nQ = 2, max.it = 20)
#'
#' # Perform prediction - use ?predict.traveltime for details
#' single_trip <- subset(tripset, tripID==2700)
#' pred <- predict.traveltime(fit, single_trip,single_trip$time[1])
#' hist(pred)
#' mean(pred)
#' sum(single_trip$traveltime)
#' 
#' ?traveltimeHMM      # for help on traveltimeHMM, the estimation function
#' ?predict.traveltime # for help on predict.traveltime, the prediction function
#' }
#' 
#' @references
#' {Woodard, D., Nogin, G., Koch, P., Racz, D., Goldszmidt, M., Horvitz, E., 2017.  Predicting travel time reliability using mobile phone GPS data.  Transportation Research Part C, 75, 30-44.}
#'
#' @importFrom stats dnorm rnorm
#' @export
traveltimeHMM <- function(logspeeds = NULL, trips = NULL, timeBins = NULL, linkIds = NULL, data = NULL,
                          nQ = 1L, model = c("HMM", "trip-HMM","trip","no-dependence"),
                          tol.err = 10, L = 10L, max.it = 20L, verbose = FALSE,
                          max.speed = NULL, seed = NULL, tmat.p = NULL,
                          init.p = NULL) {

    # SECTION A - Parameter validation and related processing
  
    # A.0 Either 'data' or all of 'logspeeds', 'trips', 'timeBins' and 'linkIds' must be provided.
    # Otherwise we stop.  Once this is checked for, we extract those fields from 'data' if required.
    frameAvailable <- !is.null(data)
    singleParamsAvailable <- !is.null(logspeeds) & !is.null(trips) & !is.null(timeBins) & !is.null(linkIds)
    
    if(!xor(frameAvailable, singleParamsAvailable))
      stop("Either 'data' or all of 'logspeeds', 'trips', 'timeBins' and 'linkIds' must be provided.")
    
    if(frameAvailable) {
      logspeeds <- data$logspeed
      trips <- data$tripID
      timeBins <- data$timeBin
      linkIds <- data$linkID
    }
    
    # param <- list(...)  # Put this line back into code later only if required.
    # Ideally all parameters should be explicit in the interface.  EG 2019/05/29
  
    # A.1 Stop if 'logspeeds', 'trips' and 'timeBins' have different sizes
    if (length(logspeeds) != length(trips) || length(logspeeds) != length(timeBins)) 
        stop("Parameter vectors for 'logspeeds', 'trips' and 'timeBins' are not equal in length!")
  
    # A.2 Parameter 'model'
    # Get model from "HMM", "trip-HMM","trip" or "no-dependence"; stop if invalid
    model <- tryCatch(match.arg(model),error=function(cond){
      stop("Parameter 'model' should be one of HMM, trip-HMM, trip, or no-dependence")
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
    
    # A.4 Parameter 'max.speed'
    # Set default value if not specified in call
    if(is.null(max.speed)){
        message('max.speed is not specified, setting at default value: 130 km/h')
        max.speed = 130
    }
    
    # We count occurrences of speeds in excess of max.speed.  If any, we then
    # report the percentage of such occurrences in a warning to the user.
    speedings = which(logspeeds > log(max.speed/3.6))
    if (length(speedings) > 0) 
        warning(paste("Many observations are higher than speed limit (130km/h)!, about", 
            paste0(round(100 * length(speedings)/length(logspeeds), 2), "% of observations."), 
            " It is advised to remove these observations."),  immediate.=TRUE, call.=FALSE)

    # SECTION B - Data preparation
    
    # Convert 'trips' and 'linkIds' to factors if required
    #if (!is.factor(trips)) 
    trips = factor(trips) 
    #if (!is.factor(linkIds)) 
    linkIds = factor(linkIds)

    nB <- length(unique(timeBins))  # time bins
    nQ2 <- nQ^2 # Compute nQ2, the square of nQ, which will be used often
    nObs <- length(logspeeds) # nObs now contains the number of observations
    nlinks <- nlevels(linkIds) # nlinks now contains the number of distinct links
    nTrips <- nlevels(trips) # nTrips now contains the number of distinct trips
    
    # timeFactor = interaction(timeBins, lex.order = TRUE) # Replaced by next line
    timeFactor = factor(timeBins)
    linkTimeFactor = interaction(linkIds, timeBins, lex.order = TRUE)  # factor of link id and time bin, using lexicographic order.
    obsId = as.numeric(linkTimeFactor) # CAUTION with that approach.  EG 2019/05/30
    tripId = as.numeric(trips) # CAUTION with that approach.  EG 2019/05/30

    # Preparation of imputation variables - refers to eq. (6) in Woodard et al.
    # We need to handle the following three cases:
    # 1. When a linkTimeFactor has fewer than L observations (counting all of them)
    # 2. When a linkTimeFactor has fewer than L initial state observations only
    # 3. When a linkTimeFactor has only initial states (so we cannot compute P(state_{k-1}, state_k | Obs}))
    # In all cases we will later impute based on the timeFactors, not taking care of the link.
    
    # Case 1
    nFactors = sapply(split(1:nObs, linkTimeFactor), length)       # Number of observations per linkTimeFactor
    linksLessMinObs = which(nFactors < L) # vector of indices of all linkTimeFactors for which we need to impute 
                                          # because of insufficient data items (<L)
    # "^[^.]+.", from start of line ^, remove everything by a dot [^.], up to the .
    indexLinksLessMinObs <- sapply(gsub("^[^.]+.", "", names(linksLessMinObs)),
                                  function(r) which(r == levels(timeFactor))) # Vector of indices to timeFactor
                                                                              # for each linkTimeFactor just identified

    # Case 2
    init_ids <- seq_along(trips)[!duplicated(trips)]  # Vector of indices of the initial state observation for each trip
    count_init <- sapply(split(1:nlevels(trips), linkTimeFactor[init_ids]),length) # vector of counts of initial
                                                                                   # state observations for each linkTimeFactor
    init_L <- which(count_init < L) # vector of indices of all linkTimeFactors with insufficient occurrence of initial occurrences
    index_init_L <- sapply(gsub("^[^.]+.", "", names(init_L)),
                           function(r) which(r == levels(timeFactor))) # vector of indices to timeFactor
                                                                       # for each linkTimeFactor just identified
                                                                       # (initial occurrence only)
  

    # Case 3
    only_init = which(sapply(split(1:nObs, linkTimeFactor), length) == count_init) # get indices of linkTimeFactors for
                                                                                   # unique trips, i.e. those for which
                                                                                   # the number of occurrences equals
                                                                                   # the number of initial occurrences
    index_only_init <- sapply(gsub("^[^.]+.", "", names(only_init)),
                             function(r) which(r == levels(timeFactor))) # Vector of indices to timeFactor
                                                                         # for each linkTimeFactor just identified...

    # Travel speed variables:
    # mu_speed and var_speed are the mean and variance matrices in eq. (3) in Woodard et al.
    # Each matrix has nQ columns and as many rows as there are linkTimeFactors.
    # At first, all rows are identical within each matrix:
    # # for mu_speed, column values are 0..nQ-1
    # # for var_speed, column values are 1..nQ
    mu_speed = matrix(1:nQ - 1, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    var_speed = matrix(1:nQ, ncol = nQ, nrow = nB * nlinks, byrow = TRUE)
    
    # Markov transition matrices and initial states per link:
    # refer to eq. (4) in Woodard et al.
    # init refers to small 'gamma' whilst tmat refers to big 'gamma'.
    
    # default initial values init.0 and tmat.0
    # WARNING: those values might be superseded, see later.  
    init.0 <- rep(1,nQ)/nQ
    tmat.0 <-rep(rep(1, nQ)/nQ, nQ)
    
    # Here we supersede the values just set for init.0 and tmat.0 if seed is passed
    if(!is.null(seed) && is.numeric(seed)){
        set.seed(seed)
        init.0 <- runif(nQ) # gives a vector of random values between 0 and 1
        init.0 <- init.0/sum(init.0) # normalizes the vector
        tmat.0 <- matrix(runif(nQ^2), nQ, nQ) # gives a matrix of random values between 0 and 1
        tmat.0 <- c(t(tmat.0/rowSums(tmat.0))) # normalizes each row and transforms into vector
    }
    # Here we attempt to supersede the value for the transition matrix if one is passed
    if(!is.null(tmat.p)){
        if(is.matrix(tmat.p) && dim(tmat.p) == c(nQ, nQ) && rowSums(tmat.p) == rep(1,nQ)){
            tmat.0 <- c(t(tmat.p))
        }else
            warning('Initial transition values are not used, tmat must be nQ x nQ with rows summing to 1', immediate.=TRUE, call.=FALSE)
    }
    
    # Here we attempt to supersede the value for the initial state vector if one is passed
    if(!is.null(init.p)){
        if(length(init.p) == nQ && sum(init.p) == 1){
            init.0 <- init.p
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
    
    # Modelling the trip effect logE: see eq. (1) from Woodard et al.
    if(grepl('trip', model)) tau2 = 1 else tau2 = 0 # set initial value for tau2: see eq. (5) from Woodard et al.
    
    # We create the vector logE of length nTrips as follows:
    # if model is "trip-HMM" or "trip" then each Ei has a random value of mean 0 and variance tau2;
    # otherwise each Ei has a value of zero.
    if(grepl('trip', model)) logE <- rnorm(nTrips, 0, sqrt(tau2)) else logE <- numeric(nTrips)
    
    ###
    # Beginning implementation of Algorithm 1 (which we call Algo1 hereafter) in Woodard et al.
    #
    ## Initialize variables - TO DOCUMENT FURTHER - EG 2019/06/03
    iter = 0  # t = 0 in step 1 of Algo1

    initNew <- tmatNew <- mu_speedNew <- var_speedNew <- probStates <- NULL
    logE_new <- logE
    tstart = Sys.time()  # starting time
    
    # Definition of error function as the absolute value of the difference between both parameters.
    # This will serve to evaluate the exit criteria from the while loop below.
    
    error.fun <- function(A, B) if (is.null(A) || is.null(B)) 0 else sum(abs(A - B))
    
    # Message to user:
    message('Model ', model, ' with ', nTrips, ' trips over ',
            nlinks, ' roads and ', nB, ' time bins...')
    
    repeat { # beginning of while loop at step 2 of Algo1
        if (grepl("HMM", model)) { # execute this block only if model is "HMM" or "trip-HMM"
          
            # probTran is an nObs x nQ^2 matrix having each row represent P(Obs_k | state_k ) * P(state_k | state_{k-1}).
            # Each row of probTran is such that matrix(probTran[i,], ncol=nQ, nrow=nQ, byrow =TRUE) is the probability
            # above in matrix format.  See documentation for 'forwardback'. 
            probTran = tmat[obsId, ] * pmax(dnorm(logspeeds, logE[tripId] + mu_speed[obsId,],
                                                  var_speed[obsId, ]), 0.001)[, rep(1:nQ, nQ)]
            
            # fb is a list of alpha and beta (respectively forward and backward probability row-wise vectors,
            # corresponding to vector 'fv' and a vector of 'b's).
            # Function 'forwardback' is called once for each trip ID with a subset of probTran
            # as first argument, and a specific value of initQ as second argument.     
            fb = tapply(1:nObs, trips, function(r)
                forwardback(probTran[r,], init[obsId[r[1]], ]))
            
            # We extract alpha and beta.
            alpha = matrix(unlist(lapply(fb, function(r) r$alpha), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # forward prob. normalized
            beta = matrix(unlist(lapply(fb, function(r) r$beta), use.names = FALSE), 
                ncol = nQ, byrow = TRUE)  # backward prob. normalized
            
            # We get the elementary product of alpha and beta (thus smoothing), then normalize by row.
            # The new matrix 'probStates' corresponds to phi_i,k(q) at step 3 of Algo1.
            probStates = normalizeR(alpha * beta + 1e-05, m = nObs, n = nQ)
            
            # We calculate the new initial state probability estimates, i.e.
            # An m x nQ matrix of probabilites for the m levels (linkTimeFactor).
            initNew = initial_est(probStates, linkTimeFactor, init_ids)
            if (length(index_init_L)) { # If there are any values to impute, then proceed
                                        # using only time factors (i.e. without information on links).
                init_impute = initial_est(probStates, timeFactor, init_ids) # impute for each timeFactor...
                initNew[init_L, ] = init_impute[index_init_L, ]             # and then replace all rows concerned
                                                                            # in the initial state probabilities matrix.
            }

            # We calculate 'probJointStates', the conditional joint transition probability P(state_k, state_{k-1} | Obs)
            # which corresponds to psi_i,k(q', q) at step 3 of Algo1.
            # We need first to remove the last row of 'probTran' and the first row of 'beta'.
            probJointStates = normalizeR(alpha[-nObs, rep(1:nQ, each = nQ)] *
                                             (probTran * beta[, rep(1:nQ, nQ)])[-1, ], m = nObs - 1, nQ2)
            probJointStates[init_ids[-1] - 1, ] <- 0  # Setting the overlapping multiplication between trips to 0

            # We compute the new transition matrix, i.e.
            # An m x nQ^2 matrix of probabilites for the m levels (linkTimeFactor).and imputing < threshold with second factor
            tmatNew = tmat_est(probJointStates, probStates, init_ids, linkTimeFactor)
            if (!is.null(indexLinksLessMinObs)) { # If there are any values to impute, then proceed
                                                  # using only time factors (i.e. without information on links).
                tmat_impute = tmat_est(probJointStates, probStates, init_ids, timeFactor) # impute for each time factor...
                tmatNew[linksLessMinObs, ] <- tmat_impute[indexLinksLessMinObs ]  # and then replace all rows concerned
                                                                                  # in the transition matrix.
            
                if (length(only_init)) # If there are trips with only one occurrence...
                    tmatNew[only_init, ] <- tmat_impute[index_only_init, ] # ... then impute using time factor.
            }
            
        }
        
        # We compute the mean and variance for the links (first two equations of step 4 in Algo1)
        meanSig = gaussian_param_by_factor(logspeeds - logE[tripId], linkTimeFactor, probStates)
        # We perform imputation on time bin only (just as before) for states with fewer than L factors.
        if (!is.null(indexLinksLessMinObs)) {
            meanSigAverage = gaussian_param_by_factor(logspeeds - logE[tripId], timeFactor, probStates)
            meanSig$mean[linksLessMinObs, ] = meanSigAverage$mean[indexLinksLessMinObs, ]
            meanSig$sigma2[linksLessMinObs, ] = meanSigAverage$sigma2[indexLinksLessMinObs, ]
        }
        
        # We enforce Woodard et al.'s restriction mu_j,b,q-1 <= mu_j,b,q (p. 34) by
        # swapping the mu_s (and sigma2_s as well) of each and every observation if required.
        # TO DO: Discuss the relevance of such an approach for enforcing the restriction.  EG 2019/06/14
        mu_speedNew = meanSig$mean
        var_speedNew = meanSig$sigma2
        # ord = order_states(meanSig$mean)
        # if (!is.null(ord$order)) {
        #     mu_speedNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) mu_speedNew[ord$toSort[r], 
        #         ord$order[r, ]], USE.NAMES = FALSE))
        #     var_speedNew[ord$toSort, ] = t(sapply(1:length(ord$toSort), function(r) var_speedNew[ord$toSort[r], 
        #         ord$order[r, ]], USE.NAMES = FALSE))
        # }
        
        # Calculating E-effect (trip specific effect parameters) --> check later by executing a trip model.  EG 2019/06/14
        if(grepl('trip', model)){
            ## Calculating E-trip variance - last equation of step 4 of Algo1
            tau2 <- .colSums(logE^2, m = nTrips, n=1)/nTrips

            ## Updating E-trip per trip
            dummy <- if(is.null(probStates)) 1/var_speedNew[obsId,] else  probStates/var_speedNew[obsId,]
            a <- .rowSums(dummy, m = nObs, n = nQ)
            h <- .rowSums(dummy * mu_speedNew[obsId,], m = nObs, n = nQ)
            dummy1 <- vapply(split(logspeeds * a - h, trips), function(r) .colSums(r, m = length(r), n=1), numeric(1), USE.NAMES = FALSE)
            dummy2 <- vapply(split(a, trips), function(r) .colSums(r, m = length(r), n=1), numeric(1), USE.NAMES = FALSE)
            logE_new <-  dummy1 / (1/tau2 + dummy2)
        }
        
        # Calculating || ThetaNew - Theta || (step 2 of Algo1)
        # In the following, A is Theta = Theta^(t-1) and B is ThetaNew = Theta^(t).
        # This will serve to compute iter_error for the purpose of exiting (or not) the loop.
        # In this version we bundle into Thetas objects of all kinds: mu_s, sigma^2_s, small_gammas, big_gammas, logE.
        # Each object has equal weight for the purpose of computing iter_error.  Isn't the result meaningless
        # from an optimization point of view.  TO DO: investigate.  EG 2019/06/14
        A = c(mu_speed[-linksLessMinObs,], var_speed[-linksLessMinObs,])
        B = c(mu_speedNew[-linksLessMinObs,], var_speedNew[-linksLessMinObs,])
        if(grepl('HMM', model)){
            A = c(A, init[-init_L,], tmat[-unique(linksLessMinObs,only_init),])
            B = c(B, initNew[-init_L,], tmatNew[-unique(linksLessMinObs,only_init),])
        }
        if(grepl('trip',model)) { A = c(A, logE); B = c(B,logE_new)}
        iter_error = error.fun(A, B)
        

        # We re-position all parameters and increment the counter.
        tmat <- tmatNew
        init <- initNew
        mu_speed <- mu_speedNew
        var_speed <- var_speedNew
        logE <- logE_new
        iter <- iter + 1

        # Announcing the expected computation length
        if(iter==1 && (!is.null(max.it) || verbose)){
            if(!is.null(max.it)){
                message('Expected completion of ', max.it, ' iterations in ',
                        format((max.it-1) * (Sys.time() - tstart), digits = 3))
            }else{
                message('Expected time per iteration ',
                        format((Sys.time() - tstart), digits = 3))
            }
        }
        
        # Announcing the "error", i.e. the distance between Theta and ThetaNew, which is
        # to decrease from one iteration to another.
        if(verbose)
            message(round(iter_error,2), " error in iteration ", iter, " @ ", format(Sys.time() - tstart, digits = 3))

                
        # Breaking loop on convergence
        if (iter_error <= tol.err) {
            message("Parameters converged at iteration " , iter-1)
            break
        }
        
        # Breaking loop on exhaution of iterations
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
                logE = logE,
                nQ = nQ,
                nB = nB,
                nObs = nObs,
                model = model)
    class(obj) <- append(class(obj),"traveltimeHMM", after=0)
    invisible(obj)
}

