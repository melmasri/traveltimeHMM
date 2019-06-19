to7daybins <-function(rules){
    if(typeof(rules[[1]])=='list' ){
        wdayrules = lapply(0:6,function(d)
            t(sapply(rules, function(r)
                if(d %in% r$days) 
                    data.frame(start = r$start, end=r$end,
                               tag = r$tag,
                               stringsAsFactors = FALSE))))
    }else{
        wdayrules = lapply(0:6,function(d)
            if(d %in% rules$days) data.frame(start = rules$start, end = rules$end,
                                             tag = rules$tag,
                                             stringsAsFactors = FALSE))
    }
    return(wdayrules)
}

#' @export
rules2timebins<-function(rules){
    time2min <- function(t){
        s = as.POSIXlt(t, format = "%H:%M")
        s$hour*60 + s$min
    }
    
    if(!is.list(rules))
        stop('rules must be a list (possibly of lists) of the format list(start,end,days,tag)!')
    if(typeof(rules[[1]])!='list'){
        if(!is.null(rules$days) & !all(rules$days %in% 0:6))
            stop('days must be in 0, 1, ..., 6, (0 = Sunday)')
        converted = rules
        converted$start = time2min(converted$start)
        converted$end = time2min(converted$end)
    }else{
        lapply(rules, function(r)
            if(!is.null(r$days) & !all(r$days %in% 0:6))
                stop('days must be in 0, 1, ..., 6, (0 = Sunday)'))
        converted = lapply(rules, function(r){
            list(start = time2min(r$start),
                 end   = time2min(r$end),
                 days = r$days,
                 tag = r$tag)
        })
    }
    ruletable =to7daybins(converted)
    return(
        time_bins_functional (
            function(t){
                ## A humanly readable function for defining rush hour only time bins
                ## returns bins of time
                day = as.POSIXlt(t)$wday
                tt = time2min(t)
                sapply(1:length(t), function(k){
                    r = ruletable[[day[k] + 1]]
                    ## weekday Mon-Fri
                    if(!is.null(unlist(r)))
                        if(length(r)==1){
                            if(tt[k] >= r['start'] &  tt[k] < r['end']) return(unlist(r['tag'][[1]]))
                        }else{
                            a =which(tt[k] >= r[,'start'] &  tt[k] < r[,'end'])
                            if(length(a)) return(unlist(r[a,'tag'][[1]]))
                        }
                    return("Other")
                })
            }
        )
    )
}

#' @export
time_bins_functional<-function(time_bin_readable_function = time_bins_readable , period = c('hours', 'minutes')){
    ## A functional function that constructs a list with bin names for each hour of
    ## the week using a humanly readable time function
    ## time 0 is Sunday 0 to 59 min AM, before 1AM
    period <- match.arg(period) # Select period, default is 'hours'
    Sun0 = as.POSIXlt("2018-08-26 00:00:00.1 EDT") # Define origin as last Sunday in August 2018 right after midnight
                                                   # Specific date is unimportant but needs to be a Sunday.

    # We define a vector of time points (in seconds) which determines time slices of length '1-hour' or '1-minute'.
    # Each time point defines the beginning of a slice.
    if(period=='hours') 
        tslice  = (0:(24*7-1))*3600 # case 'hours' : 1-hour time slices
    else 
        tslice  = (0:(24*60*7-1))*60 # case 'minutes' : 1-minute time slices
    
    # Get a vector of time bins for each time point thus defined
    time.bins.per.slice = sapply(tslice, function(r) time_bin_readable_function(Sun0+r))
    
    # Return a function that gives the time bin of a given date using appropriate slicing.
    if(period == 'hours'){
        return(
            function(t){
                # time 0 is Sunday 0 to 59 min AM, before 1AM
                t = as.POSIXlt(t)
                time.bins.per.slice[t$hour + 1 + t$wday*24]
            }
        )
    }else{
        return(
            function(t){
                # time 0 is Sunday 0 to 59 min AM, before 1AM
                t = as.POSIXlt(t)
                time.bins.per.slice[t$hour + 1 + t$wday*24 + t$min]
            }
        )
    }
}

time_bins_readable <- function(t){
    # A humanly readable function that returns bins of time
    day = as.POSIXlt(t)$wday # Get day of week: from Sunday (0) to Saturday (6)
    h = as.POSIXlt(t)$hour # Get hour
    sapply(1:length(t), function(k){ # Return time bin corresponding to hour and day of week
        if(day[k] %in%  1:5){
            ## weekday Mon-Fri
            if(h[k] >= 7 & h[k] < 9) return('MorningRush')
            if(h[k] >= 15 & h[k]< 18) return('EveningRush')
            if((day[k]==5 & (h[k]>= 20 | h[k]<6)) | (day[k] %in% 1:4 & (h[k]>=19 | h[k]<6)))
                return ('EveningNight')
            return('Weekday')
        }else{
            ## weekend Sunday Saturday
            if(day[k]==6 & (h[k] >= 21 | h[k] < 9)) return("EveningNight")
            if(day[k]==0 & (h[k] >= 19 | h[k] < 9)) return("EveningNight")
            return("Weekendday")
        }
    })
}

time_bins <- time_bins_functional()
