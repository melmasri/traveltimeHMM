#' \code{cleanGPS} cleans raw GPS data
#'

#' @details Cleans raw GPS data by removing trips with outside the indicated parameters.
cleanGPS <- function(raw,
                     min.gps.obs = 50,
                     max.idle.time = 4,
                     max.time.between.gps =  120,
                     min.speed  = 5,
                     median.speed = 5,
                     min.max.speed = 9,
                     min.distance = 1000,
                     min.endpoint.speed = 10
                     ){
    ## ##################################################
    ## ### Settings 
    ## ##################################################
    ## min.gps.obs = 50                         # in units
    ## max.idle.time = 4                # mins
    ## median.speed  = 5                             # m/s
    ## min.max.speed = 9                             # m/s
    ## min.distance = 1000                   # in meters
    ## min.speed=5                                   # in km/h
    ## max.time.between.gps = 120                    # sec
    ## min.endpoint.speed  = 10                     # km/h
    ## ##################################################
    
    ## #-------------------------------------------------- mesc functions
    ## Needed functions and libraries
    require(data.table)
    require(geosphere)
    require(mclust)
    options(digits.secs=3)
    dist_geo<-function(x,y, min){
        z= cbind(x,y)
        z = rbind(min,z)
        distCosine(z)
    }
    cumsum_time<-function(t) t -min(t)
    rush_hour<-function(t) (t>= 6 & t <=9) | (t>=15 & t <=18)
    ## start processes time
    

    ## # -------------------------------------------------- sanity checks
    ## checking if data.{frame, table}
    if(!is.data.frame(raw) || !is.data.table(raw))
        stop('raw must be a data.frame or data.table object')

    ## convert to data.table
    if(!is.data.table(raw)){
        trips = data.table(raw)
    }else{
        trips = raw
    }

    ## checking column names
    col.names<-c('trip', 'timestamp', 'speed', 'lat', 'long', 'roadId')
    if(all(colnames(trips) %in% col.names))
        stop('raw column names must be (trip, timestamp, speed, lat, long, roadId)')
    
    ## #-------------------------------------------------- Descriptive statistics
    sTime = Sys.time()
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))
    print('trips per hour')
    trips[,.(hour=as.numeric(format(time[1],format='%H'))),by=trip][,.(no.trips=length(unique(trip))),by=hour][order(hour)]
    
    
##################################################
    ## removing trips with small number of observations
    print(paste('removing trips with less than', min.gps.obs, 'observations'))
    print(paste('total to be removed', sum(trips[,.N <= min.gps.obs, by =trip]$V1)))
    trips = trips[trips[,.I[.N > min.gps.obs], by =trip]$V1]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
    ## adjusting data
    trips[, geoDist:= dist_geo(long, lat,c(long[1], lat[1])), by=trip]
    trips[, location:=cumsum(geoDist), by=trip]
    trips[, speedEst:= 3.6*c(0,diff(location)/as.numeric(diff(time))), by=trip]
    trips[, speed:= 3.6*speed]


###############################################
### decomposing trips
    print('Sum of consecutive stopping time')
    groups_<-function(r){
        if(length(r)==1) return(1)
        aux = which(diff(r)>1)
        if(length(aux)==0)  return(rep(1,length(r)))
        n = length(r)
        ind = rep(length(aux)+1, n)
        for(i in aux)
            ind[1:i] = ind[1:i]-1
        return(ind)
    }

    aux = trips[geoDist<(min.speed/3.6)][,.(time, g=groups_(roadId),speed, speedEst,roadId),by=trip]
    aux = aux[,.(tdiff = difftime(time[.N], time[1], units='secs'),id1 = roadId[1], id2=roadId[.N]),by=c('trip', 'g')]

    print('percetage of trips with no move for the following durations')
    n=nrow(aux)
    t = seq(0,10,0.5)
    rbind(t, no = colSums(aux[, outer(tdiff, t*60, '<=')])/n)

    print(paste('Decomposing trips with no move (<',min.speed/3.6,'m/s, or ',min.speed,'k/m)for more than', max.idle.time, 'mins in between'))
    aux = aux[tdiff > max.idle.time*60]
    max.trip.id = 3000000
    print(paste('number of trips to be decomposed', nrow(aux), 'trips'))
    print(paste('new trips start with id', max.trip.id))

    ls = list()
    for(i in 1:nrow(aux)){
        a = trips[, which(trip==aux[i,trip])]
        id1 = a[which(trips[a,id==aux[i,id1]])]
        id2 = a[which(trips[a,id==aux[i,id2]])]
        max.trip.id = max.trip.id + 1
        ls[[i+1]]<-list( I = id2:a[length(a)], id = max.trip.id)
    }
    ## updating trip ids
    aux = sapply(ls, function(r) cbind(r$I, rep(r$id, length(r$I))))
    aux = do.call('rbind', aux)
    t = trips$trip
    t[aux[,1]]<-aux[,2]
    trips[, trip:=t]

    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

### removing 0 geoDist from top and bottom of each trip
    trim.zero.ends<-function(g, s=min.speed){
        a = which(g>s)
        if(length(a)==0) return(FALSE)
        if(length(a)==1) return(FALSE)
        if(length(a)>1) return(max((min(a)-1),1):max(a))
    }

### removing speed < min.speed from top and bottom of each trip
    trips= trips[,  .SD[trim.zero.ends(speed, min.endpoint.speed-1e-5)],by=trip]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))
    trips= trips[order(trip,time)]

### removing 0 geoDist from top and bottom of each trip
    trips= trips[,  .SD[trim.zero.ends(geoDist,0)],by=trip]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))
    trips= trips[order(trip,time)]

##################################################
    ## removing trips with small number of observations after decomposition
    print(paste('removing trips with less than', min.gps.obs, 'observations'))
    print(paste('total to be removed', sum(trips[,.N <= min.gps.obs, by =trip]$V1)))
    trips = trips[trips[,.I[.N > min.gps.obs], by =trip]$V1]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
### breaking trips with time diff larger than max.time.between.gps
    print(paste('breaking trips with obs time difference more than', max.time.between.gps/60, 'mins'))
    groups.per.time.diff<-function(r){
        aux = which(diff(r)>= max.time.between.gps)
        if(length(aux)==0)  return(rep(1,length(r)))
        n = length(r)
        ind = rep(length(aux)+1, n)
        for(i in aux)
            ind[1:i] = ind[1:i]-1
        return(ind)
    }

    aux = trips[ , .(g = groups.per.time.diff(time), time, roadId) , by = trip]
    aux = aux[, .SD[max(g)>1], by=trip]
    aux = aux[, .(id1= roadId[1], id2=roadId[.N], g), by=c('trip', 'g')]
    print(paste('Decomposing trips with no move for more than',  max.time.between.gps/60, 'mins in between'))
    print(paste('number of trips to be decomposed', aux[,1,by=trip][,sum(V1)], 'trips'))

    max.trip.id = 4000000
    print(paste('new trips start with id', max.trip.id))
    ls = list()
    for(i in 1:nrow(aux)){
        a = trips[, which(trip==aux[i,trip])]
        id1 = a[which(trips[a,id==aux[i,id1]])]
        id2 = a[which(trips[a,id==aux[i,id2]])]
        max.trip.id = max.trip.id + 1
        ls[[i+1]]<-list( I = id2:a[length(a)], id = max.trip.id)
    }

    ## updating trip ids
    aux = sapply(ls, function(r) cbind(r$I, rep(r$id, length(r$I))))
    aux = do.call('rbind', aux)
    t = trips$trip
    t[aux[,1]]<-aux[,2]
    trips[, trip:=t]

    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
    ## removing trips with small number of observations after decomposition
    print(paste('removing trips with less than', min.gps.obs, 'observations'))
    print(paste('total to be removed', sum(trips[,.N <= min.gps.obs, by =trip]$V1)))
    trips = trips[trips[,.I[.N > min.gps.obs], by =trip]$V1]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
    ## adjusting data can calculating other quntities
    trips[, geoDist:= dist_geo(long, lat,c(long[1], lat[1])), by=trip]
    trips[, location:=cumsum(geoDist), by=trip]
    trips[, speedEst:= 3.6*c(0,diff(location)/as.numeric(diff(time))), by=trip]

### removing speed less than min.endpoint.speed 
    trips= trips[,  .SD[trim.zero.ends(speed,min.endpoint.speed-1e-5)],by=trip]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))
    trips= trips[order(trip,time)]

##################################################
    ## adjusting data can calculating other quntities
    trips[, geoDist:= dist_geo(long, lat,c(long[1], lat[1])), by=trip]
    trips[, location:=cumsum(geoDist), by=trip]
    trips[, speedEst:= 3.6*c(0,diff(location)/as.numeric(diff(time))), by=trip]

### removing 0 geoDist from top and bottom of each trip
    trips= trips[,  .SD[trim.zero.ends(geoDist,0)],by=trip]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))
    trips= trips[order(trip,time)]

##################################################
    ## removing trips with small number of observations after decomposition
    print(paste('removing trips with less than', min.gps.obs, 'observations'))
    print(paste('total to be removed', sum(trips[,.N <= min.gps.obs, by =trip]$V1)))
    trips = trips[trips[,.I[.N > min.gps.obs], by =trip]$V1]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
### removing trips with median speed and min.max.speed
    print(paste('Removing trips with median estimated speed <', median.speed, 'm/s or max estimated speed <', min.max.speed, 'm/s.'))
    aux = trips[speedEst>0, median(speedEst) >= 3.6*median.speed & max(speed) >= 3.6*min.max.speed, by =trip]
    print(paste('Total removed', sum(!aux$V1), 'trips'))
    trips = trips[trip %in% aux$trip[aux$V1]]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
    ## removing trips short driving distance
    ## preliminary cleaning
    print(paste('removing trips with driving distance less than', min.distance, 'meters.' ))
    print(paste('total to be removed', sum(trips[,max(location) - min(location)<= min.distance, by=trip]$V1)))
    trips = trips[trips[,.I[max(location) - min(location)> min.distance], by=trip]$V1]
    print(paste('No. of trips:', nrow(trips[, .N, by=trip])))

##################################################
    ## removing walkers
    cluster_based_walkers <-function(loc, speed, min.speed){
        a = sum(diff(loc)[speed[-1]<min.speed])
        a / (max(loc) - min(loc))
    }

    v = trips[, cluster_based_walkers(location, speedEst, min.speed) ,by=trip]
    cluster = Mclust(v$V1, G = 20)
    cluster$parameters$mean
    class = v$trip[which(cluster$classification %in% which(cluster$parameters$mean>0.30))]

    print(paste('removing walker/runer trips with 30% of the time the speed less than', min.speed, 'km/h'))
    print(paste('Total to be removed', length(class), 'trip(s)'))
    trips = trips[!(trip %in% class)]
    cat('No. of trips:', nrow(trips[, .N, by=trip]), '\n')


##################################################
### writing file
    print('writing file')
    trips[, location:=NULL]
    trips[, geoDist:=NULL]
    trips[, speedEst:=NULL]
    trips = trips[, .(roadId, trip, timestamp, speed, long, lat)]
    message('Total time ',difftime(Sys.time(), sTime, units='mins'), ' mins.')
    return(invisible(trips))
}
