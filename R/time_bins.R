

to7daybins <-function(rules){
    if(typeof(rules[[1]])=='list' ){
        wdayrules = lapply(0:6,function(d)
           lapply(rules, function(r)
               if(d %in% r$days)
                   list(start = r$start, end=r$end,
                        tag = r$tag, stringsAsFactors = FALSE)
                  ))
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
    time_bins = time_bins_functional (
            function(t){
                ## A humanly readable function for defining rush hour only time bins
                ## returns bins of time
                day = as.POSIXlt(t)$wday
                tt = time2min(t)
                getTags <-function(x){
                    if(!is.null(x))
                        if(tt[k] >= x['start'] & tt[k] < x['end']) x['tag']
                }
                sapply(1:length(t), function(k){
                    r = ruletable[[day[k] + 1]]
                    ## weekday Mon-Fri
                    ## If the whole list is null
                    if(is.null(unlist(r))) return ('Other')
                    ## otherwise search along lists
                    tag = unlist(lapply(r, getTags), use.names=FALSE)
                    if(!is.null(tag)) tag[1] else  return("Other")
                })
            }
        )
    
    return(time_bins)
}

#' @export
time_bins_functional<-function(time_bin_readable_function = time_bins_readable , period = c('hours', 'minutes')){
    ## A functional function that constructs a list with bin names for each hour of
    ## the week using a humanly readable time function
    ## time 0 is Sunday 0 to 59 min AM, before 1AM
    period <- match.arg(period)
    Sun0 = as.POSIXlt("2018-08-26 00:00:00.1 EDT")
    if(grepl('hour', period))
        tslice  = (0:(24*7-1))*3600
    if(grepl('min', period))
        tslice  = (0:(24*60*7-1))*60
    time.bins.per.slice = sapply(tslice, function(r) time_bin_readable_function(Sun0+r))
    if(period == 'hours'){
        return(
            function(t){
                ## time 0 is Sunday 0 to 59 min AM, before 1AM
                t = as.POSIXlt(t)
                time.bins.per.slice[t$hour + 1 + t$wday*24]
            }
        )
    }else{
        return(
            function(t){
                ## time 0 is Sunday 0 to 59 min AM, before 1AM
                t = as.POSIXlt(t)
                time.bins.per.slice[t$hour + 1 + t$wday*24 + t$min]
            }
        )
    }
}

time_bins_readable <- function(t){
    ## A humanly readable function that 
    ## returns bins of time
    day = weekdays(t)
    h = as.POSIXlt(t)$hour
    sapply(1:length(t), function(k){
        if(day[k] %in%  c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")){
            ## weekday Mon-Fri
            if(h[k] >= 7 & h[k] < 9) return('MorningRush')
            if(h[k] >= 15 & h[k]< 18) return('EveningRush')
            if((day[k]=='Friday'& (h[k]>= 20 | h[k]<6)) | (day[k] %in% c("Monday", "Tuesday", "Wednesday", "Thursday") & (h[k]>=19 | h[k]<6)))
                return ('EveningNight')
            return('Weekday')
        }else{
            ## weekend Sunday Saturday
            if(day[k]=="Saturday" & (h[k] >= 21 | h[k] < 9)) return("EveningNight")
            if(day[k]=="Sunday" & (h[k] >= 19 | h[k] < 9)) return("EveningNight")
            return("Weekendday")
        }
    })
}

time_bins <- time_bins_functional()
