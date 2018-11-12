rm(list=ls())
library(data.table)
load('../ready2model_simple.RData')

devtools::load_all()

analyze.prediction<-function(pred, obs, file = NULL,plot=TRUE, ...){
    dt = merge(pred, obs, by='trip')
    pointEst = dt[, .(
        pred = mean(predTT),
        obs = obsTT[1],
        timeBin = timeBins[1],
        len = len[1],
        ARE =  abs(mean(predTT) -obsTT[1])/obsTT[1],
        AR = abs(mean(predTT) - obsTT[1]),
        bias =log(mean(predTT)) - log(obsTT[1] ),
        Covered = sign(prod(quantile(predTT, probs = c(0.025, 0.975)) - obsTT[1]))==-1,
        QuantileRange = diff(quantile(predTT, probs = c(0.025,0.975)))
    ), by=trip]

    errors = pointEst[, .(
        geoARE = exp(mean(log(ARE))),
        meanAE = mean(AR),
        bias = mean(bias),
        meanARE = mean(ARE),
        coverege = mean(Covered),
        covergeInterval = mean(QuantileRange),
        coverageOfObs = mean(QuantileRange/obs)
    ),]
    errorsTimeBin = pointEst[, .(
        geoARE = exp(mean(log(ARE))),
        meanAE = mean(AR),
        bias = mean(bias),
        meanARE = mean(ARE),
        coverege = mean(Covered),
        covergeInterval = mean(QuantileRange),
        coverageOfObs = mean(QuantileRange/obs)
    ),by = timeBin]
    print('Errors over all')
    print(errors)
    print('Errors over timeBins')
    print(errorsTimeBin)

    coverage_intervals = seq(0,1, 0.05)
    names(coverage_intervals)<-paste0('alpha_', coverage_intervals)
    quant = dt[,{
        pred = mean(predTT);
        obs = obsTT[1];
        coverage = lapply(coverage_intervals, function(a) diff(quantile(predTT, probs = c(a/2, 1-a/2), names=FALSE)));
        covered =  lapply(coverage_intervals, function(a) sign(prod(quantile(predTT, probs = c(a/2, 1-a/2)) - obsTT[1]))==-1);
        list(pred = pred, obs=obs, coverage = list(coverage), covered = list(covered))
    },by=trip]
    
    return(invisible(list(errors = errors,
                          quantiles = quant,
                          empirical.coverage = rowMeans(sapply(quant$covered, unlist))
                        , interval.width = rowMeans(sapply(quant$coverage, unlist) ))))
}

##################################################
### Simulation run 
## data loading
aux = which(tt.trip.link$logspeed > log(130/3.6))
if(length(aux)>0){
    warning(paste('Many observations are higher than speed limit! 130km/h', 'removing ', length(aux), 'observation, %', round(100*length(aux)/nrow(tt.trip.link),2)))
    tt.trip.link= tt.trip.link[-aux]
}

seed = 12346
set.seed(seed)
no.test = 1500
test.trips = sample(unique(tt.trip.link$trip), no.test)
newlinkIds = tt.trip.link[(trip %in% test.trips)][, linkidrel]
aux = setdiff(newlinkIds, tt.trip.link[!(trip %in% test.trips), linkidrel])
length(aux)
## remove trips with single links 
test.trips = setdiff(test.trips, tt.trip.link[linkidrel %in% aux][, unique(trip)])
## write.csv(test.trips, file = 'test_trips.csv', row.names = FALSE)
length(test.trips)

tt.trip.link.train = tt.trip.link[!(trip %in% test.trips)]
## Estimate model parameters
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel, maxItra=45, model ='HMM', nQ = 3)
## predict for held-out-set
pred = tt.trip.link[trip %in% test.trips, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
## 
obs  = tt.trip.link[trip %in% test.trips, .(obsTT = sum(tt), timeBins=timeBins[1], len = sum(length)), by=trip]

analyze.prediction(pred, obs)

## saveRDS(aux$empirical.coverage, file='woodard_hmm_empirical-coverge.RDS')
## saveRDS(aux$interval.width, file='woodard_hmm_interval-width.RDS')

###################################################
### No - dependence
seed = 12346
set.seed(seed)
test.trips = unlist(read.csv(file = 'test_trips.csv'))

tt.trip.link.train = tt.trip.link[!(trip %in% test.trips)]
## Estimate model parameters
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel, maxItra=45, model ='no-dependence')
## predict for held-out-set

pred = tt.trip.link[trip %in% test.trips[1:100], .(predTT= predict.traveltime(est, linkidrel, length,time[1], n=100)), by=trip]
##
obs  = tt.trip.link[trip %in% test.trips[1:400], .(obsTT = sum(tt)), by=trip]

aux = analyze.prediction(pred, obs)

saveRDS(aux$empirical.coverage, file='woodard_no-dependece_empirical-coverge.RDS')
saveRDS(aux$interval.width, file='woodard_no-dependece_interval-width.RDS')




## plottig
fullmodelEmprical =readRDS(file='woodard_hmm_empirical-coverge.RDS')
fullmodelInterval = readRDS(file='woodard_hmm_interval-width.RDS')


nodepEmprical = readRDS(file='woodard_no-dependece_empirical-coverge.RDS')
nodepInterval = readRDS(file='woodard_no-dependece_interval-width.RDS')

plot(fullmodelEmprical, nodepEmprical)

cbind(fullmodelEmprical, nodepEmprical)

pdf('interval.pdf')
plot(1:length(nodepInterval), nodepInterval, col='blue', type='b')
lines(1:length(nodepInterval), fullmodelInterval, col='red', type='b')
dev.off()


##################################################
#### in sample test

### markov model 
seed = 12346
set.seed(seed)
test.trips = unlist(read.csv(file = 'test_trips.csv'))

## Estimate model parameters
est = traveltimeHMM(tt.trip.link$logspeed, tt.trip.link$trip, tt.trip.link$timeBins, tt.trip.link$linkidrel, maxItra=45, model ='HMM', nQ = 2)
## predict for held-out-set
pred = tt.trip.link[trip %in% test.trips[1:10], .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
## 
obs  = tt.trip.link[trip %in% test.trips[1:10], .(obsTT = sum(tt)), by=trip]
aux = analyze.prediction(pred, obs)

### no-dep
seed = 12346
set.seed(seed)
test.trips = unlist(read.csv(file = 'test_trips.csv'))
## Estimate model parameters
est = traveltimeHMM(tt.trip.link$logspeed, tt.trip.link$trip, tt.trip.link$timeBins, tt.trip.link$linkidrel, maxItra=45, model ='no-dependence')
## predict for held-out-set
pred = tt.trip.link[trip %in% test.trips, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
## 
obs  = tt.trip.link[trip %in% test.trips, .(obsTT = sum(tt)), by=trip]
aux = analyze.prediction(pred, obs)



## plotting travel time per links
load('woodard_pred.RData')
library(data.table)
plot.tt.links<-function(lk,timeBin, file=NULL, ...){
    sampling.size = 1000
    if(!is.null(file))
        bitmap(file, height = 6, width=6, units = 'in', res=300)
    par(mfrow = rep(ceiling(length(lk)/2) , 2))
    for(t in lk){
        par(mar = c(5,5,1,1) + 0.1)
        N = tt.trip.link[timeBins == timeBin & linkidrel == t,.N]
        d = links[linkidrel == t, length]
        osm = links[linkidrel ==t, osm]
        
        tt.trip.link[timeBins == timeBin & linkidrel == t,
                     hist(tt, breaks = 50, freq=FALSE,
                          xlab = 'travel time (s)', cex.lab=1.5, lty=1, lwd=2,
                          main = paste('osm', osm ,
                              ',sample', N, ',', round(d,0),'m'))]
        a = which(levels(est$factors) == paste(t, timeBin, sep='.'))
        q = sample.int(est$nQ,sampling.size, replace = TRUE, prob =  est$init[a,])
        s = rnorm(sampling.size, est$mean[a,q], est$sd[a, q])
        lines(density(d * exp(-s)), lwd = 2)
    }
    if(!is.null(file)) dev.off()
}


## Top links 100 links
set.seed(12345)
timeBin = 'EveningRush'
top.links = sample(tt.trip.link[timeBins == timeBin, .N,  by = linkidrel][order(N, decreasing=TRUE)][1:100, linkidrel], 4)
sampling.size = 1000
plot.tt.links(top.links, 'Top100.png')


## lower links
set.seed(2344)
timeBin = 'EveningRush'
bottom.links = sample(tt.trip.link[timeBins == timeBin, .N,  by = linkidrel][order(N)][N>10][1:3000, linkidrel], 4)
sampling.size = 1000
plot.tt.links(bottom.links, 'bottom100.png')

## travel time for a path
set.seed(12345)
timeBin = 'EveningNight'
tt.trip.link[timeBins == timeBin, .N,  by = linkidrel][order(N, decreasing=TRUE)]

## trip 1
trip1= c(157685071,25794864 , 39276632, 39276633)
linkids = c(32086,6510,20674,20676,20675)
maxlen = links[linkidrel %in% linkids, sum(length)]
trips1 =   tt.trip.link[timeBins == timeBins & linkidrel %in% linkids, .(N=length(linkids), len=sum(length)), by=trip][N==length(linkids) & len == maxlen, trip]


len = links[linkidrel %in% linkids, length][match(linkids, links[linkidrel %in% linkids,linkidrel])]
      
tt.trip.link[timeBins == timeBins & linkidrel %in% linkids & trip %in% trips1, unique(length), by=linkidrel]
tt = tt.trip.link[timeBins == timeBin & linkidrel %in% linkids & trip %in% trips1, sum(tt), by=trip]$V1
tt.week = tt.trip.link[timeBins == 'Weekday' & linkidrel %in% linkids & trip %in% trips1, sum(tt), by=trip]$V1
plot.tt.links(linkids, timeBin, file = 'Trip_travelTime_2km_per_link.png')
plot.tt.links(linkids, 'Weekday', file = 'Trip_travelTime_2km_per_link_weekday_trip_sigle.png')

png('Trip_travelTime_2km.png', height = 6, width=6, units = 'in', res=300)
par(mar = c(5,5,1,1) + 0.1)
hist(tt, breaks = 60, freq=FALSE,
     xlab = 'travel time (s)', cex.lab=1.5, lty=1, lwd=2,
     main = paste('sample', length(trips1), ',', round(maxlen,0),'m'), col=rgb(1,0,0,0))
hist(tt.week, add=TRUE, col=rgb(0,0,1,0.5), freq=FALSE, breaks = 60, cex.lab=1.5, lty=1, lwd=2)

stime = tt.trip.link[timeBins==timeBin, min(time)]
lines(density(predict.traveltime.HMM(est, linkids,len, starttime = stime,n =5000)), lwd = 5,  col='red')

stime= tt.trip.link[timeBins=='Weekday'][1, time]
lines(density(predict.traveltime.HMM(est, linkids,len, starttime = stime,n =5000)), lwd = 5, col='blue')
dev.off()






