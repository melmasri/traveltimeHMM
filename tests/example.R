rm(list=ls())
devtools::document()
devtools::load_all()

## load data
data(trips)
q1trips = trips[1:20000,]
rules = list(
    list(start='6:30',  end= '9:00',  days = 1:5, tag='MR'),
    list(start='15:00', end= '18:00', days = 2:5, tag='ER')
)

time_bins <- rules2timebins(rules)


data(trips)
fit <- traveltimeHMM(trips$logspeed, trips$trip, trips$timeBin, trips$linkId, nQ = 2, max.it = 5)
single_trip <- subset(trips, trip==2700)
pred <- predict.traveltime(fit, single_trip, time_bins.fun = time_bins)
hist(pred, main ='travel time distribution for trip 2700', freq =FALSE)      # histogram of prediction samples
abline(v = mean(pred), lty=2, lwd=2)
mean(pred)      # travel time point estimate
abline(v=sum(single_trip$traveltime), lty=2, lwd=2, col='red')    # observed travel time

est = traveltimeHMM(trips$logspeed, trips$trip, trips$timeBin, trips$linkId, nQ = 2, verbose=TRUE)


single_trip$time[1]


## Cleaning
devtools::load_all()
trips = readRDS('../traveltime/dataTableDatabase/Quebec_trips.rds')
library(data.table)
trips[, 1, trip][, sum(V1)]

trips[, timestamp :=time]
trips[, roadId:=id]
trips[, ':='(time=NULL, id=NULL)]

n = nrow(trips)
cleaned = tapply(1:n, trips$trip, function(r){
    print(trips[r[1], trip])
    cleanGPS_trip(trips[r,], trips[r[1], trip])
})



table(unlist(sapply(cleaned, function(r)
    lapply(r, function(s)
        if(s$tag=='dirty') s$comment else  'clean'))))
    

cleanedsubset = lapply(cleaned, function(r){
    aux = which(sapply(r, function(s) s$tag) == "clean")
    if(length(aux)){
        r = r[aux]
        k = lapply(r, function(s) {s$trip[, tripId:=s$trip_id];s$trip})
        rbindlist(k)
        return(k)
    }
})

cleanedsubset  = cleanedsubset[-which(sapply(cleanedsubset, is.null))]

cleanedsubset = lapply(cleanedsubset, function(r){
    rbindlist(r)
})

cleaned = rbindlist(cleanedsubset, use.names=TRUE)

cleaned[, 1, tripId][,sum(V1)]
cleaned[, trip:=tripId]
cleaned[, tripId:=NULL]

saveRDS(cleaned, 'trips_clean.rds')

### Creating full dataset
devtools::load_all()
rules = list(
    list(start='6:30', end= '9:00', days = 1:5, tag='MR'),
    list(start='15:00', end= '18:00', days = 1:5, tag='ER')
)
time_bins <- rules2timebins(rules)

### Input settings
fileName = 'trips_osm.rds'
fileMeta = 'trips_osm_meta.rds'
fileGeo  = 'osm_geometry.rds'

### testing covariance assumption

devtools::load_all()
trips = readRDS('trips.rds')
library(data.table)

z = 1:15
u = qchisq(.975, df=14)
l = qchisq(1-.975, df=14)
largen = trips[ , .N, trip][N>20, trip]
X = trips[trip %in% largen,
      {
      y = abs(acf(exp(logspeed), plot=FALSE, lag.max = 15)[[1]][-1]);
      x= predict(lm(y~1 + z))
      sign(sum((y-x)^2/y) - u)
  }, by = trip]

X[, mean(V1==-1)]

trips[trip==23052, acf(exp(logspeed), plot=FALSE, lag.max =15)]


