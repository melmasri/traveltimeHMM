
time_bins_readable <- function(t){
    ## A humanly readable function that 
    ## returns bins of time
    day = weekdays(t)
    h = hour(t)
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

time_bins_functional<-function(time_bin_readable_function){
    ## A functional function that constructs a list with bin names for each hour of
    ## the week using a humanly readable time function
    ## time 0 is Sunday 0 to 59 min AM, before 1AM
    Sun0 = as.POSIXlt("2018-08-26 00:00:00.1 EDT")
    weekhours = 0:(24*7-1)
    time.bins.per.week.hour = sapply(weekhours, function(r) time_bin_readable_function(Sun0+r*3600))
    function(t){
        ## time 0 is Sunday 0 to 59 min AM, before 1AM
        whour = hour(t)+ 1 + (wday(t)-1)*24
        time.bins.per.week.hour[whour]
    }
}

time_bins <- time_bins_functional(time_bins_readable)
