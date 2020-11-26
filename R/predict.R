#' Predict the travel time for a trip using a \code{traveltimeHMM} model object
#' 
#' \code{predict.traveltimeHMM} performs a point prediction by simulation using parameter estimates provided by a \code{traveltimeHMM} model object.
#' Prediction can be performed for a single trip only.
#' 
#' The function begins by validating and, if required, replacing the value of the parameter \code{logE}
#' (see explanation alongside \code{logE} in the \emph{Arguments} section).  It then transfers execution
#' to the appropriate function according to the selected model: \code{predict.traveltimeHMM} for
#' models of the \code{HMM} family, or \code{predict.traveltimeHMM.no_dependence} otherwise.
#'  
#' @param object A model object (a list) provided through the execution of function \code{timetravelHMM}.
#'   The list includes information on model as well as estimates for its parameters.
#'   See \code{timetravelHMM} man page.
#' @param tripdata A data frame of road links with information on each
#'   link's traversal.  Columns minimally include objects 'linkID' and 'length',
#'   and the latter must have the same length.  Rows must be in chronological order.
#'   The program assumes that the sequence of road links forms a coherent and feasible
#'   path.  \emph{No verification is performed to that effect}.
#' @param starttime The start date and time for the very first link of the trip,
#'   in POSIXct format.  Default is the current date and time.
#' @param n Number of samples.  Default is 1000.
#' @param logE Point estimate of trip effects.  \code{logE} normally needs to be a numerical vector of size \code{n}.
#'   If a single numerical value is supplied, it will be replicated into a vector.  If \code{logE} is \code{NULL}
#'   the function will use either a vector of simulated values (if the model is from the \code{trip} family),
#'   or a vector of \code{0} otherwise.  Default is \code{NULL}.  NOTE: when simulating values for the
#'   vector, the value for \eqn{\tau} is taken from the model object.
#' @param time_bins.fun A functional to map real time to specified time bins, see `?rules2timebins`.
#' @param ... not used.
#'
#' @return \code{predict.traveltimeHMM} returns a numerical vector of size \code{n} representing the point prediction of total travel time, in seconds, for each run.
#'
#' @examples
#' \dontrun{
#' data(tripset)
#' 
#' # Fit a model - use ?traveltimeHMM for details
#' fit <- traveltimeHMM(tripset$logspeed, tripset$tripID,
#'                      tripset$timeBin, tripset$linkID, nQ = 2, max.it = 10)
#' 
#' # Perform a prediction for trip #2700 using the fitted model.
#' single_trip <- subset(tripset, tripID==2700)
#' 
#' # We need to supply the time stamp of the very first link traversal (third parameter)
#' pred <- predict(fit, single_trip,single_trip$time[1])
#'
#' hist(pred)                     # histogram of prediction samples
#' mean(pred)                     # travel time point estimate
#' sum(single_trip$traveltime)    # observed travel time
#' 
#' ?traveltimeHMM      # for help on traveltimeHMM, the estimation function
#' ?predict.traveltimeHMM # for help on predict.traveltimeHMM, the prediction function
#' }
#' @references
#' {Woodard, D., Nogin, G., Koch, P., Racz, D., Goldszmidt, M., Horvitz, E., 2017.  Predicting travel time reliability using mobile phone GPS data.  Transportation Research Part C, 75, 30-44.}
#' @export
predict.traveltimeHMM<-function(object, tripdata, starttime = Sys.time(),  n = 1000, logE = NULL,time_bins.fun = time_bins, ... ){
  
    # We first perform basic checks.  'tripdata' must be a list, data frame or data table
    # that minimally includes objects 'linkID' and 'length', the latter having
    # the same length.
    if(!is.list(tripdata))
      stop('tripdata must be a list, data.frame or data.table')
    if(!all(c('linkID', 'length') %in% names(tripdata)))
      stop('tripdata must have objects named linkID and length, corresponding order of travelled links and their lengths')
    if(length(tripdata$linkID)!=length(tripdata$length))
      stop('length of objects do not match!')
  
    # Models of the HMM family ('HMM', 'trip-HMM') are handled by function 'predict.traveltime.HMM'
    # whilst others are handled by function 'predict.traveltimeHMM.no_dependence'
    # (both functions are below).
    if(grepl('HMM', object$model))
      predict.traveltimeHMM.HMM(object, tripdata, starttime, n, logE, time_bins = time_bins.fun, ...)
    else
      predict.traveltimeHMM.no_dependence(object, tripdata , starttime, n, logE, time_bins = time_bins.fun, ...)
}

#' Predict the travel time for a trip using a \code{traveltimeHMM} model object that is not of the HMM family
#' @keywords internal
#' 
#' \code{predict.traveltimeHMM.no_dependence} performs a point prediction by simulation using parameter estimates provided by a \code{traveltimeHMM} model object that is not of the \code{HMM} family (see man page for \code{predict.traveltimeHMM}).
#' 
#' The function implements Algorithm 2 from Woodard et al., 2017.  However, the state transition matrix
#' and initial state probability vector are not handled as they were not generated at the estimation stage.
#'  
#' @param object A model object (a list) provided through the execution of
#' function \code{timetravelHMM} for a \code{trip} or \code{no-dependence} model type.
#' The list includes information on model as well as estimates
#' for its parameters.  See \code{timetravelHMM} man page.
#' @param tripdata A data frame of road links with information on each
#'   link's traversal.  Columns minimally includes objects 'linkID' and 'length',
#'   and the latter must have the same length.  Rows must be in chronological order.
#' @param starttime The start date and time for the very first link of the trip, in POSIXct format.
#' @param n Number of samples. Default is 1000.
#' @param logE Point estimate of trip effects, in the form of a numerical vector of size \code{n}.
#' @param time_bins a functional map between real time and time bins, see `?rules2timebins`.
#' @param ... not used.
#' 
#' @return \code{predict.traveltimeHMM.no_dependence} returns a vector of size \code{n} of representing the point prediction of total travel time, in seconds, for each run.
#' 
#' @importFrom stats rnorm runif
#' @references
#' {Woodard, D., Nogin, G., Koch, P., Racz, D., Goldszmidt, M., Horvitz, E., 2017.  Predicting travel time reliability using mobile phone GPS data.  Transportation Research Part C, 75, 30-44.}
#' @export
predict.traveltimeHMM.no_dependence <- function(object, tripdata, starttime, n = 1000, logE = NULL, time_bins = time_bins, ...) {
    linkIds = tripdata$linkID # Contains IDs of all links for a given trip
    len = tripdata$length # Contains the length (in km) of each link in 'linkIds'
    logE <- getValidE(object, logE, n) # Get a valid vector for 'logE'; see comments in function for details.

    # Get link+time factor ID for link on top of list and start time supplied
    fact = paste(linkIds[1], time_bins(starttime), sep = ".")
    id = which(levels(object$factors) == fact) # Find ID of factor corresponding to the string just created

    ###
    # Beginning implementation of Algorithm 2 (which we call Algo2 hereafter) in Woodard et al.
    #
    ## Initialize variables
    
    # Generate a vector of size n of random speeds with mu and sigma supplied by 'object'.
    # This corresponds to step 9 of Algo2, but only for the first iteration.
    speed = rnorm(n, object$mean[id, ], object$sd[id, ])
    tt = len[1] * exp(-speed - logE) # n-sized vector of total travel time
                                     # on first link: length * speed (adjusted for trip effect)
    for (k in 2:length(linkIds)) { # Loop for each link after the first
        # Get n-sized vector of link+timebin factors
        fact = as.factor(paste(linkIds[k], time_bins(starttime + tt), sep = "."))

        # Get vector of link+time factor IDs for link supplied and resulting time
        id = sapply(levels(fact), function(s) which(levels(object$factors) == s, useNames = FALSE), USE.NAMES = FALSE)
        ind = as.numeric(fact) # Get n-sized vector of factor IDs relative to vector 'id'

        # Step 9 of Algo2: generate a vector of size n of random speeds
        # with mu and sigma supplied by 'object'.
        speed = rnorm(n, object$mean[id[ind], ], object$sd[id[ind], ])

        # Step 10 of Algo2: generate a vector of size n of traversal times for the current link
        tt = tt + len[k] * exp(-speed - logE) # total travel time for the current link
                                              # = length * speed (adjusted for trip effect)
    }
    tt
}


#' Predict the travel time for a trip using a \code{traveltimeHMM} model object of the HMM family
#' @keywords internal
#' 
#' \code{predict.traveltime.HMM} performs a point prediction by simulation using parameter estimates provided by a \code{traveltimeHMM} model object of the \code{HMM} family (see man page for \code{predict.traveltimeHMM}).
#' 
#' The function implements Algorithm 2 from Woodard et al., 2017, including its handling
#' of the state transition matrix and initial state probability vector.
#'  
#' @param object A model object (a list) provided through the execution of
#' function \code{timetravelHMM} for a \code{trip-HMM} or \code{HMM} model type.
#' The list includes information on model as well as estimates
#' for its parameters.  See \code{timetravelHMM} man page.
#' @param tripdata A data frame of road links with information on each
#'   link's traversal.  Columns minimally includes objects 'linkID' and 'length',
#'   and the latter must have the same length.  Rows must be in chronological order.
#' @param starttime The start date and time for the very first link of the trip, in POSIXct format.
#' @param n Number of samples. Default is 1000.
#' @param logE Point estimate of trip effects, in the form of a numerical vector of size \code{n}.
#' @param time_bins a functional map between real time and time bins, see `?rules2timebins`.
#' @param ... not used.
#' 
#' @return \code{predict.traveltimeHMM.HMM} returns a vector of size \code{n} of representing the point prediction of total travel time, in seconds, for each run.
#' 
#' @importFrom stats rnorm runif
#' @references
#' {Woodard, D., Nogin, G., Koch, P., Racz, D., Goldszmidt, M., Horvitz, E., 2017.  Predicting travel time reliability using mobile phone GPS data.  Transportation Research Part C, 75, 30-44.}
#' @export
predict.traveltimeHMM.HMM <- function(object, tripdata, starttime, n, logE, time_bins = time_bins, ...) {
    linkIds = tripdata$linkID # Contains IDs of all links for a given trip
    len = tripdata$length # Contains the length (in km) of each link in 'linkIds'
    logE <- getValidE(object, logE, n) # Get a valid vector for 'logE'; see comments in function for details.
    nQ = object$nQ # Get number of states from object.

    # Get link+time factor ID for link on top of list and start time supplied
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
                         sep = "."), useNames = FALSE)
    if(length(id)==0) # Stop if no such factor is found
      stop("No link and time combination corresponds to those supplied at the beginning.  Stopping.")

    ###
    # Beginning implementation of Algorithm 2 (which we call Algo2 hereafter) in Woodard et al.
    #
    ## Initialize variables
    
    # Step 3 of Algo2: generate a vector of states 1:nQ of size n.
    # The probability of each state is fixed and provided in object$init
    # for the first link and start time supplied).
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])

    # Generate a vector of size n of random speeds with mu and sigma supplied by 'object'.
    # This corresponds to step 9 of Algo2, but only for the first iteration.
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    
    tt = len[1] * exp(-speed - logE) # n-sized vector of total travel time
                                     # on first link: length * speed (adjusted for trip effect)
    if (length(linkIds) > 1)
        for (k in 2:length(linkIds)) { # Loop for each link after the first
          
            # Get n-sized vector of link+timebin factors
            # (This saves about 50ms for 117 routes (.4 ms per route))
            fact = as.factor(paste(linkIds[k], time_bins(starttime + tt), sep = "."))

            # Get vector of link+time factor IDs for link supplied and resulting time
            id = sapply(levels(fact), function(s) which(levels(object$factors) == 
                s, useNames = FALSE), USE.NAMES = FALSE)
            if(length(id)==0) # Stop if no such factor is found
              stop(paste("No link and time combination corresponds to those supplied at k =", k, ".  Stopping."))
            
            ind = as.numeric(fact) # Get n-sized vector of factor IDs relative to vector 'id'
            indIds = (1:length(id) - 1) * nQ
            ## indIds is suppsed to be the first index of each unique factor, if length(id) =1 , indIds = c(0,0), if length(id)=2 then inIds = c(0,2)
            ## the logic behind indIds, is that if nQ =2, then the first two rows for tmat2 corresponds to the first unique factor, the second 2 rows to the second unique factor
            ## etc, .. the number of unique factors is length(id)
            ## creating tmat takes about 10ms for 2 rows or the full matrix to speed up create
            ## tmat2 before the loop or post-estimation for the whole Quebec dataset, for 1294
            ## obs it takes about an extra 11sec if at the top of the loop
            
            # When there's more than 1 id,
            # aren't we mixing data from different factors???  EG 2019/06/19
            tmat2 = matrix(c(t(object$tmat[id, ])), ncol = nQ, byrow = TRUE)
            ## tmat2 is an nQ x length(id)
            ## the first two rows are the transition of levels(fact)[1], the seond 2 rows are the transitions of levels(fact)[2], and so on.
            # if nQ =2 , lenth(id)=1 this is 2x2 matrix tmat[1, ] = (q1,q1), (q1,q2), (q2,q1), (q2,q2),
            ## tmat2  = (q1,q1) , (q2, q1), 2nd row (q1,q2), (q2,q2), tmat2 has rows equal to unique factors of fact
            tmat2 = t(apply(tmat2, 1, cumsum))  # here tmat2 = c(total prob. to go to q1, total prob. to go to q2) = c(p(q1), p(q2)) if length(id)=1 an nQ=2
            tmat2 = tmat2[indIds[ind] + Qk, ]   # QK = 1, 2, 1, 2, 2, 1, 2, ... is the old states
            ## since indIds is the index of the top row of the transition matrix for each factor
            ## then, asumme ind correspond to the 2nd unique factor
            ## then inIds[ind] would return 2, meaning row 3, and 4 of tmat2 is transition matrix of this fact
            ## if the old state Qk=1, then 2+1, returns row 3, if old state is 2, this returns row 4.
            ## Another example, assume that ind = 4th unique factor, then indIds[4] = 3*2 = 6, meaning the 4th factor has transition matrix in levels 7 and 8.
            ## if the old state is 1 then  tmat2[indIds[ind] + Qk, ] returne row 7 = indIds[4] + 1, otherwise returns row 8
            ## and so on.
            ## QK below samples a new state based on the new tmat2, which correspods now to 1 row per unqiue factor, and it is the row of the transition matrix of the old state.
            ## now that I think about it, I am not sure about this line tmat2 = t(apply(tmat2, 1, cumsum)), need to be checked.            Qk = max.col(runif(n) < tmat2, "first")

            # Step 9 of Algo2: generate a vector of size n of random speeds
            # with mu and sigma supplied by 'object'.
            speed = rnorm(n, object$mean[cbind(id[ind], Qk)], object$sd[cbind(id[ind], Qk)])

            # Step 10 of Algo2: generate a vector of size n of traversal times for the current link
            tt = tt + len[k] * exp(-speed - logE) # total travel time for the current link
                                        # = length * speed (adjusted for trip effect)
        }
    tt
}

