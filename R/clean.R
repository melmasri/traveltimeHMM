#' \code{cleanGPS} cleans raw GPS data
#'
#'
#' @details Cleans raw GPS data by removing trips with outside the indicated parameters.
#' @export
cleanGPS_trip<-function(trip,trip_id= NULL,verbose= FALSE, ...){
    listed_trip = cleanGPS_trip_internal(trip,trip_id,verbose, ...)
    ## to flatten the list
    ## listed_trip comes as a tree and the unlist_trip makes it just a one level tree
    ## Flatten it
    unlist_trip<-function(a){
        new_list  = list()
        unlist_<-function(L){
            if(length(L)==0) return(NULL)
            if(!is.null(L$tag)){
                new_list[[length(new_list)+1]]<<-L
                return(NULL)
            }
            if(!is.null(L[[1]]$tag)){
                new_list[[length(new_list)+1]]<<-L[[1]]
                L[[1]]<-NULL
                unlist_(L)
            }else{
                if(length(L)>1){
                    list(unlist_(L[[1]]),unlist_(L[-1]))
                }else{
                    unlist_(L[[1]])
                }
            }
        }
        x = unlist_(a)
        return(new_list)
    }
    flat = unlist_trip(listed_trip)
    return(flat)
}

cleanGPS_trip_internal<-function(trip,trip_id= NULL,verbose= FALSE, ...){
    ## ##################################################
    ## ### Input parameters
    ## ##################################################
    ## min.gps.obs = 50                 # in units
    ## max.idle.time = 4                # mins
    ## median.speed  = 20               # km/h
    ## min.max.speed = 35               # km/h
    ## min.distance = 1000              # meters
    ## min.speed = 5                    # km/h
    ## max.time.between.gps = 120       # sec
    ## min.endpoint.speed  = 10         # km/h
    ## ##################################################
    ## #-------------------------------------------------- flow of cleaning
    ## 1. order trips by trip and timestamp
    ## 2. calculate estimated speed and consecutive location distance
    ## 3. trimming endpoints of trips such that the trips starts and ends when speed crosses min.endpoint.speed
    ## 4. removing trips if less than median speed, or min.max speed
    ## 5. decompose trips with idle time larger than max.idle.time, return to 1.
    ## 6. decompose trips with time difference larger than max.time.between.gps, return to 1.
    ## 7. return if clean
    
    arg = list(...)
    min.gps.obs =  if(is.null(arg$min.gps.obs)) 50 else arg$min.gps.obs
    max.idle.time = if(is.null(arg$max.idle.time)) 4 else arg$max.idle.time
    max.time.between.gps = if(is.null(arg$max.time.between.gps)) 120 else arg$max.time.between.gps
    min.speed  = if(is.null(arg$min.speed)) 5 else arg$min.speed
    median.speed = if(is.null(arg$median.speed)) 20 else arg$median.speed
    minmax.speed = if(is.null(arg$minmax.speed)) 35 else arg$minmax.speed
    min.distance = if(is.null(arg$min.distance)) 1000 else arg$min.distance
    min.endpoint.speed =if(is.null(arg$min.endpoint.speed)) 10 else arg$min.endpoint.speed
        
    comment = NULL
    tag = NULL
    ## #-------------------------------------------------- required functions
    require(geosphere)
    require(data.table)
    dist_geo<-function(x,y, min){
        z = rbind(min,cbind(x,y))
        distHaversine(z)
    }
    cumsum_time<-function(t) t -min(t)
    ## # -------------------------------------------------- sanity checks
    ## checking if data.{frame, table}
    if(!is.list(trip))
        stop('trip must be a list.')
    if(!is.data.table(trip))
        trip = data.table(trip)
    ## checking column names

    names.var<-c('timestamp', 'speed', 'lat', 'long', 'roadId')
    if(!all(names.var %in% names(trip)))
        stop('trip names must be (timestamp, speed, lat, long, roadId)')
    
    if(length(trip$roadId)<=1){
        comment = 'a single observation trip, please remove!'
        if(verbose)
            warning(comment)
        return(list(trip, tag = 'dirty', comment = comment))
    }
    
    ## #--------------------------------------------------
    ## adjusting data
    trip = trip[order(timestamp)]
    trip[, geoDist := dist_geo(long, lat, c(long[1], lat[1]))]
    trip[, location := cumsum(geoDist)]
    trip[, speedEst := 3.6*c(0,diff(location)/as.numeric(diff(timestamp)))] # in k/h

    ## #--------------------------------------------------
    ## trimming endpoints of a trip for min.speed
    trim.zero.ends<-function(g, s=min.speed){
        a = which(g>s)
        if(length(a)==0) return(FALSE)
        if(length(a)==1) return(FALSE)
        if(length(a)>1) return(max(min(a)-1,1):max(a)) # keeping the first obs time
    }
    ## removing speed < min.speed from top and bottom of each trip
    trip = trip[,  .SD[trim.zero.ends(speedEst, min.endpoint.speed-1e-5)]]
    trip= trip[order(timestamp)]

    ## return if trip has one observatoin
    if(length(trip$roadId)<=1){
        comment = 'a single observation trip, please remove!'
        if(verbose)   warning(comment)
        trip[, ':='(geoDist=NULL, speedEst=NULL, location=NULL)]
        return(list(trip = trip,trip_id = trip_id, tag = 'dirty', comment = comment))
    }


    ## #--------------------------------------------------
    ## decomposing trips based on max.idle.time
    groups_<-function(r){
        if(length(r)==1) return(1)
        aux = which(diff(r)>1)
        if(length(aux)==0)  return(rep(1,length(r)))
        n = length(r)
        ind = rep(length(aux)+1, n)
        for(i in aux)                   # this part needs a speed-up
            ind[1:i] = ind[1:i]-1
        return(ind)
    }

    if(!is.null(min.speed) && !is.null(max.idle.time)){
        trip[,reluId:=.I]                   # reluId is unique per row
        aux =trip[geoDist<(min.speed/3.6)][,.(timestamp, g=groups_(reluId),
            speedEst,roadId, reluId)]
        if(nrow(aux)>1){
            aux = aux[,.(tdiff = difftime(timestamp[.N],timestamp[1],units='secs'),
                id1 = reluId[1], id2=reluId[.N]),by=g]
            aux = aux[tdiff > max.idle.time*60]
            if(nrow(aux)>=1){
                if(verbose){
                    message('Decomposing trip with no move (< ',min.speed,'k/m) for more than ', max.idle.time, 'mins in between')
                    message('trip decomposed into ', nrow(aux)+1, ' trips')
                }
                a = aux[, c(0,id2)]
                if(a[length(a)]!=nrow(trip)) a = c(a,nrow(trip))
                a = rep(1:(length(a)-1), diff(a))
                trip[, reluId:=NULL]
                g = a
                tag = 'D'
                return(tapply(1:nrow(trip), g, function(r)
                    cleanGPS_trip(trip[r], paste0(trip_id, tag, g[r[1]]),verbose, ...)
                              )
                       )
            }
        }
        trip[, reluId:=NULL]            # removing unwanted variable
    }
    ## breaking trips with time diff larger than max.time.between.gps
    if(!is.null(max.time.between.gps)){
        if(verbose)
            message('breaking trips with obs time difference more than ', max.time.between.gps/60, 'mins')
        groups.per.time.diff<-function(r){
            aux = which(diff(r)>= max.time.between.gps)
            if(length(aux)==0)  return(rep(1,length(r)))
            n = length(r)
            ind = rep(length(aux)+1, n)
            for(i in aux)               # this needs speed-up
                ind[1:i] = ind[1:i]-1
            return(ind)
        }
        g = groups.per.time.diff(trip$timestamp)
        if(max(g)>1){
            tag = 'D'
            return(
                return(tapply(1:nrow(trip), g, function(r)
                    cleanGPS_trip(trip[r], paste0(trip_id, tag,g[r[1]]),verbose, ...)
                              )
                       )
            )
            if(verbose)
                message('number of trips to be decomposed ',nrow(aux) , ' trips.')
        }
    }
    ## #-------------------------------------------------- removing based on median speed
    if(!is.null(median.speed) && !is.null(minmax.speed)){
        if(verbose)
            message('Removing trips with median estimated speed <', median.speed, 'km/h or max speed <', minmax.speed, 'km/h.')
        
        aux = trip[speedEst>0, median(speedEst) >= median.speed & max(speedEst) >= minmax.speed]
        if(!aux){
            comment = 'Median or max speed are low!'
            if(verbose) warning(comment)
            trip[, ':='(geoDist=NULL, speedEst=NULL, location=NULL)]
            return(list(trip = trip, trip_id = trip_id,tag='dirty', comment = comment))
        }
    }
    ## #-------------------------------------------------- min.distance
    ## removing trips short driving distance
    if(!is.null(min.distance)){
        ## preliminary cleaning
        if(verbose)
            message('removing if driving distance less than ', min.distance, 'm.' )
        aux  =trip[, max(location) - min(location)> min.distance]
        if(!aux){
            comment = 'trip shorter then min distance!'
            if(verbose) warning(comment)
            trip[, ':='(geoDist=NULL, speedEst=NULL, location=NULL)]
            return(list(trip = trip, trip_id = trip_id,tag='dirty', comment = comment))
        }
    }
    trip[, ':='(geoDist=NULL, speedEst=NULL, location=NULL)]
    ## #-------------------------------------------------- clean stage
    return(list(trip = trip, trip_id = trip_id, tag='clean'))
}

## trips[, timestamp :=time]
## trips[, roadId:=id]
## trips[, ':='(time=NULL, id=NULL)]
## trip= trips[trip==2690]
## a=  cleanGPS_trip(trip, trip_id = 2690, verbose= FALSE)

