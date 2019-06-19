
#' @export

# Function validateE returns the appropriate value for E.
# E can be passed either as a separate parameter with default value NULL, 
# or as part of 'object'.  Valid values of parameter E are considered first and override
# any value for E in 'object'.  However, any value passed (either by parameter or
# through 'object) need to be a vector of length 'n'; otherwise it is discarded
# If no valid value is found, then we either assign a vector of random values following
# a normal distribution with mean 0 and standard deviation tau (if trip model), or a
# vector of zeros (otherwise).

getValidE <- function(object, E, n) {
  # Part 1 of validation.  We check for the existence of parameter E.
  if(is.null(E)) { # If not found, we assign value in 'object' if valid.
    if('E' %in% names(object) && !is.null(object$E)) {
      if(is.numeric(object$E) && length(object$E)==n) 
        E <- object$E
      else
        message("Values E in object need to be a vector of length n.  Values for E are discarded.")
    }
  } else if(!is.numeric(E) || length(E)!=n) { # If a parameter exists, use it only if valid,
                                              # otherwise set to NULL.
    message("Values in parameter E need to be a vector of length n.  Values for E are discarded.")
    E <- NULL
  }
  # At this point, E has the appropriate value among valid parameters passed,
  # and NULL if no such value is found.
  
  # Part 2 of validation: if E is null...
  
  if(is.null(E)) {
    if(grepl('trip', object$model)) 
      E = rnorm(n, mean = 0, sd = object$tau) # ... we define a vector of random values if trip model...
    else E = 0 # ... otherwise we set to zero.
  }
  E # We return the result.
}


predict.traveltime<-function(object, data, starttime = Sys.time(),  n = 1000, E = NULL){
  
    # We first perform basic checks.  'data' must be a list, data frame or data table
    # that minimally includes objects 'linkId' and 'length', the latter having
    # the same length.
    if(!is.list(data))
      stop('data must be a list, data.frame or data.table')
    if(!all(c('linkId', 'length') %in% names(data)))
      stop('data must have objects named linkId and length, corresponding order of travelled links and their lengths')
    if(length(data$linkId)!=length(data$length))
      stop('length of objects do not match!')
  
    # Models of the HMM family ('HMM', 'trip-HMM') are handled by function 'predict.traveltime.HMM'
    # whilst others are handled by function 'predict.traveltime.no_dependence'
    # (both functions are below).
    if(grepl('HMM', object$model))
      predict.traveltime.HMM(object, data, starttime, n, E)
    else
      predict.traveltime.no_dependence(object, data , starttime, n, E)
}

predict.traveltime.no_dependence <- function(object, data, starttime, n = 1000, E = NULL) {
    linkIds = data$linkId
    len = data$length
    ## sampling E (trip-effect)
    E <- getValidE(object, E, n)
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


predict.traveltime.HMM <- function(object, data, starttime, n = 1000, E = NULL) {
    ## sampling E (trip-effect)
    linkIds = data$linkId # Contains IDs of all links for a given trip
    len = data$length # Contains the length (in km) of each link in 'linkIds'
    E <- getValidE(object, E, n) # Get a valid vector for 'E'; see comments in function for details.
    nQ = object$nQ # Get number of states from object.

    # Get link+time factor ID for link on top of list and start time supplied
    id = which(levels(object$factors) == paste(linkIds[1], time_bins(starttime), 
                         sep = "."), useNames = FALSE)
    if(length(id)==0) # Stop if no such factor is found
      stop("No link and time combination corresponds to those supplied at the beginning.  Stopping.")

    ###
    # Beginning implementation of Algorithm 2 (which we call Algo2 hereafter) in Woodard et al.
    #
    ## Initialize variables - TO DOCUMENT FURTHER - ÉG 2019/06/19
    
    # Step 3 of Algo2: generate a vector of states 1:nQ of size n.
    # The probability of each state is fixed and provided in object$init
    # for the first link and start time supplied).
    Qk = sample.int(nQ, n, replace = TRUE, prob = object$init[id, ])

    # Step 2 of Algo2: generate a vector of size n of random speeds with mu and sigma supplied by 'object'.
    speed = rnorm(n, object$mean[id, Qk], object$sd[id, Qk])
    
    tt = len[1] * exp(-speed - E) # n-sized vector of total travel time
                                  # on first link: length * speed (adjusted for trip effect)
    if (length(linkIds) > 1) 
        for (k in 2:length(linkIds)) {
          
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
            
            ## creating tmat takes about 10ms for 2 rows or the full matrix to speed up create
            ## tmat2 before the loop or post-estimation for the whole Quebec dataset, for 1294
            ## obs it takes about an extra 11sec if at the top of the loop
            
            # Clarification is needed regarding the workings of
            # the creation of tmat2 below.  When there's more than 1 id,
            # aren't we mixing data from different factors?  ÉG 2019/06/19
            tmat2 = matrix(c(t(object$tmat[id, ])), ncol = nQ, byrow = TRUE)
            tmat2 = t(apply(tmat2, 1, cumsum))
            tmat2 = tmat2[indIds[ind] + Qk, ]
            Qk = max.col(runif(n) < tmat2, "first")
            speed = rnorm(n, object$mean[cbind(id[ind], Qk)], object$sd[cbind(id[ind], Qk)])
            tt = tt + len[k] * exp(-speed - E)
        }
    tt
}

