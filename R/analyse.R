
analyze.prediction<-function(pred, obs, file = NULL,plot=TRUE, ...){
    dt = merge(pred, obs, by='trip')
    pointEst = dt[, .(pred = mean(predTT),
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
