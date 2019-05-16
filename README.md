# traveltimeHMM
Travel time prediction from GPS points using an HMM.

# Installation
```
install.packages("devtools")
devtools::install_github("melmasri/traveltimeHMM", auth_token = 'get an Auth_token') # this is for privacy measure
##devtools::install_github("melmasri/traveltimeHMM")
```
To get an auth_token see [this](https://help.github.com/en/articles/creating-a-personal-access-token-for-the-command-line).

# for an example
```R
library(traveltimeHMM)
data(trips)
?traveltimeHMM  # for help
fit <- traveltimeHMM(trips$logspeed, trips$trip, trips$timeBin, trips$linkId, nQ = 2, max.it = 20)
single_trip <- subset(trips, trip==2700)
pred <- predict.traveltime(fit, single_trip,single_trip$time[1])
hist(pred)      # histogram of prediction samples
mean(pred)      # travel time point estimate
sum(single_trip$traveltime)    # observed travel time
```
