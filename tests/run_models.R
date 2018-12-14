rm(list=ls())
library(data.table)
devtools::load_all()
load('ready2model.RData')

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
    coverage_intervals = seq(0,1, 0.01)
    names(coverage_intervals)<-paste0('alpha_', coverage_intervals)
    quant = dt[,{coverage = lapply(coverage_intervals, function(a) diff(quantile(predTT, probs = c(a/2, 1-a/2), names=FALSE)));
                 covered =  lapply(coverage_intervals, function(a) sign(prod(quantile(predTT, probs = c(a/2, 1-a/2)) - obsTT[1]))==-1);
                 list( coverage = list(coverage), covered = list(covered))
             },by=trip]
    
    pointEst = merge(pointEst, quant, by = 'trip')

    errors = pointEst[, .(N=.N,
        geoARE = exp(mean(log(ARE))),
        meanAE = mean(AR),
        bias = mean(bias),
        meanARE = mean(ARE),
        coverege = mean(Covered),
        covergeInterval = mean(QuantileRange),
        coverageOfObs = mean(QuantileRange/obs),
        interval95 =  rowMeans(sapply(coverage, unlist))[which(rowMeans(sapply(covered, unlist)) <= .95)[1]]),]

    errorsTimeBin = pointEst[, .(
        N=.N,
        geoARE = exp(mean(log(ARE))),
        meanAE = mean(AR),
        bias = mean(bias),
        meanARE = mean(ARE),
        coverege = mean(Covered),
        covergeInterval = mean(QuantileRange),
        coverageOfObs = mean(QuantileRange/obs),
        interval95 =  rowMeans(sapply(coverage, unlist))[which(rowMeans(sapply(covered, unlist)) <= .95)[1]]
    ),by = timeBin]


    coverage_intervals = seq(0,1, 0.05)
    names(coverage_intervals)<-paste0('alpha_', coverage_intervals)
    quant = dt[,{
        pred = mean(predTT);
        obs = obsTT[1];
        coverage = lapply(coverage_intervals, function(a) diff(quantile(predTT, probs = c(a/2, 1-a/2), names=FALSE)));
        covered =  lapply(coverage_intervals, function(a) sign(prod(quantile(predTT, probs = c(a/2, 1-a/2)) - obsTT[1]))==-1);
        list(pred = pred, obs=obs, coverage = list(coverage), covered = list(covered))
    },by=trip]
    
    print('Errors over all')
    print(errors)
    print('Errors over timeBins')
    print(errorsTimeBin[order(timeBin)])
    return(invisible(list(errors = errors,
                          quantiles = quant,
                          empirical.coverage = rowMeans(sapply(quant$covered, unlist))
                        , interval.width = rowMeans(sapply(quant$coverage, unlist) ))))
}

badObs = which(tt.trip.link$logspeed > log(140/3.6))
badtrips = NULL
if(length(badObs)>0){
    warning(paste('Many observations are higher than speed limit! 140km/h', 'removing ', length(badObs), 'observation, %', round(100*length(badObs)/nrow(tt.trip.link),2)))
    badtrips = tt.trip.link[badObs, unique(trip)]
    tt.trip.link= tt.trip.link[-badObs]
}

est = traveltimeHMM(tt.trip.link$logspeed, tt.trip.link$trip, tt.trip.link$timeBins, tt.trip.link$linkidrel,  model ='no-dep', max.it=200, tol.err = 10, nQ =1)
paramid = as.numeric(est$factors)
devT.test = tt.trip.link[, .(test = abs((logspeed - est$mean[paramid]))/est$sd[paramid], trip, linkidrel )][, .I[test>4]]
if(length(devT.test)){
    warning(paste('remove outliers in the lower end, approximitly', round(100*length(devT.test)/nrow(tt.trip.link),2), '% will be removed'))
    badtrips = c(badtrips,  tt.trip.link[devT.test, unique(trip)])
    badtrips = unique(badtrips)
    tt.trip.link= tt.trip.link[-devT.test]
}

## -------------------------------------------------- setting up test sets
seed = 12346
set.seed(seed)
no.test = 1500
testtrips = sample(setdiff(unique(tt.trip.link$trip), badtrips), no.test)
newlinkIds = tt.trip.link[(trip %in% testtrips)][, linkidrel]
aux = setdiff(newlinkIds, tt.trip.link[!(trip %in% testtrips), linkidrel])
## remove trips with single links 
testtrips = setdiff(testtrips, tt.trip.link[linkidrel %in% aux][, unique(trip)])
testtrips = sort(setdiff(testtrips, tt.trip.link[trip %in% testtrips][, .N==1, by = trip][V1==TRUE, trip]))
## write.csv(test.trips, file = 'test_trips.csv', row.names = FALSE)
length(testtrips)
tt.trip.link.train = tt.trip.link[!(trip %in% testtrips)][order(trip, time)]

## -------------------------------------------------- trip-HMM model
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel,  model ='trip-HMM', max.it=1000, tol.err = 10, nQ =2)
saveRDS(est, file = 'est_trip-HMM.rds')
## predict for held-out-set
pred = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
obs  = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(obsTT = sum(tt), timeBins=timeBins[1], len = sum(length)), by=trip]
aux = analyze.prediction(pred, obs)
write.csv(aux$empirical.coverage, file='woodard_trip-hmm_empirical-coverge.csv')
write.csv(aux$interval.width, file='woodard_trip-hmm-_interval-width.csv')

## -------------------------------------------------- HMM model
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel,  model ='HMM', max.it=100, tol.err = 10, nQ =2)
saveRDS(est, file = 'est_HMM.rds')
## predict for held-out-set
pred = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
obs  = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(obsTT = sum(tt), timeBins=timeBins[1], len = sum(length)), by=trip]
aux = analyze.prediction(pred, obs)
write.csv(aux$empirical.coverage, file='woodard_hmm_empirical-coverge.csv')
write.csv(aux$interval.width, file='woodard_hmm_interval-width.csv')

## -------------------------------------------------- trip model
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel,  model ='trip', max.it=1000, tol.err = 10, nQ =1)
saveRDS(est, file = 'est_trip.rds')
## predict for held-out-set
pred = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
obs  = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(obsTT = sum(tt), timeBins=timeBins[1], len = sum(length)), by=trip]
aux = analyze.prediction(pred, obs)
write.csv(aux$empirical.coverage, file='woodard_trip_empirical-coverge.csv')
write.csv(aux$interval.width, file='woodard_trip_interval-width.csv')


## -------------------------------------------------- no-dependece model
est = traveltimeHMM(tt.trip.link.train$logspeed, tt.trip.link.train$trip, tt.trip.link.train$timeBins, tt.trip.link.train$linkidrel,  model ='no-dep', max.it=100, tol.err = 10, nQ =1)
saveRDS(est, file = 'est_no-dep.rds')
## predict for held-out-set
pred = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(predTT= predict.traveltime(est, linkidrel, length,time[1])), by=trip]
obs  = tt.trip.link[trip %in% testtrips][order(trip,time)][, .(obsTT = sum(tt), timeBins=timeBins[1], len = sum(length)), by=trip]
aux = analyze.prediction(pred, obs)
write.csv(aux$empirical.coverage, file='woodard_no-dependece_empirical-coverge.csv')
write.csv(aux$interval.width, file='woodard_no-dependece_interval-width.csv')


##-------------------------------------------------- plotting empirical coverage
files = dir(pattern = 'empirical')
files = files[order(files)]
png('Coverage_per_model.png', height = 6, width=6, units = 'in', res=300)
par(mar = c(5,5,1,1) + 0.1)
aux = read.csv(files[1])[,2]
plot(seq(0, 1, by = 0.05), aux[length(aux):1],cex.lab=1.5, lty=1, lwd=2,
     xlab = 'Theoretical Coverage of Predictive Interval', ylab = 'Empirical Coverage of Test Data', type='b', ylim =c(0,1), xlim = c(0,1))
abline(a = 0, b = 1, col='black', lwd=3)
for(k in 2:length(files)){
    aux = read.csv(files[k])[,2]
    lines(seq(0, 1, by = 0.05), aux[length(aux):1],
          cex.lab=1.5, lty=1, lwd=2, col=k, type='b')
}
legend("topleft", legend = gsub('(woodard_|_empirical-coverge.csv)', '', files),
       lwd = 2, col = c(1:length(files)), cex=1.5)
dev.off()

## -------------------------------------------------- plotting interval width
files = dir(pattern = 'interval')
files = files[order(files)]
png('Interval_width_per_model.png', height = 6, width=6, units = 'in', res=300)
par(mar = c(5,5,1,1) + 0.1)
aux = read.csv(files[1])[,2]
plot(seq(0, 0.95, by = 0.05), aux[length(aux):2],cex.lab=1.5, lty=1, lwd=2,
     xlab = 'Theoretical Coverage of Predictive Interval', ylab = 'Empirical Interval Width', type='b', ylim = c(0, 900))
for(k in 2:length(files)){
    aux = read.csv(files[k])[,2]
    lines(seq(0, 0.95, by = 0.05), aux[length(aux):2],
          cex.lab=1.5, lty=1, lwd=2, col=k, type='b')
}
legend("topleft", legend = gsub('(woodard_|_interval-width.csv)', '', files),
       lwd = 2, col = c(1:length(files)), cex=1.5)
dev.off()



## -------------------------------------------------- linear model

print('linear model')
fit = lm(log(tt) ~ log(len) + timeBins, data = tt.trip.link.train[, .(tt = sum(tt), len = sum(length), timeBins = timeBins[1]), by =trip])

tt.trip.link.test = tt.trip.link[trip %in% testtrips]
obs = tt.trip.link.test[, .(obsTT = sum(tt),
    len = sum(length),
    timeBins = timeBins[1]), by = trip]
obs = obs[order(trip, timeBins)]

pred = cbind(trip = obs$trip, predTT = exp(predict(fit, newdata = obs, prediction.interval=TRUE)))
pred = data.table(pred)
analyze.prediction(pred, obs)

pred = cbind(trip = obs$trip, exp(predict(fit, newdata = obs,interval = 'pred')))
dt = data.table(merge(pred, obs, by = 'trip' ))

dt[, Covered := lwr<= obsTT & upr>=obsTT]
dt[, QuantileRange := upr - lwr]

dt[, .(coverege = mean(Covered),
       covergeInterval = mean(QuantileRange),
       coverageOfObs = mean(QuantileRange/obsTT))]

dt[, .(coverege = mean(Covered),
       covergeInterval = mean(QuantileRange),
       coverageOfObs = mean(QuantileRange/obsTT)), by = timeBins]


