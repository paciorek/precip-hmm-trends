library(lubridate)

startYr <- 1920 # Must be >= 1880
endYr <- 2021   # Must be <= 2018 (although no data past 2018-03-15)

cat("Loading ", data_fn, ".\n")
# orig <- read.csv(file.path('data','PA_station.csv'))
orig <- read.csv(file.path(data_dir, paste0(data_fn, '.csv')))


if(names(orig)[1] == 'time') { # old csv format
    names(orig) <- c('date','raw_year','season_year','month','prcp')
    endYr <- 2017
}

orig <- orig[orig$date != '1900-2-29', ]  # this day didn't actually exist (1900 not a leap year)

season_months <- list(c(12,1,2), 3:5, 6:8, 9:11, c(11,12,1,2,3,4))
season_mid <- c('1-15','4-15','7-15','10-15', '1-31')

if(season == 5) {  # need to move November to next year
    orig$season_year[orig$month %in% c(11)] <- orig$season_year[orig$month %in% c(11)] + 1
    lastNov <- which(orig$season_year == endYr+1)
    orig$prcp[lastNov] <- NA
    orig$raw_year[lastNov] <- min(orig$raw_year)
    orig$season_year[lastNov] <- min(orig$season_year)
    orig <- rbind(orig[lastNov,], orig[-lastNov,])
    substring(orig$date[orig$month %in% c(11)], 1, 4) <- as.character(orig$raw_year[orig$month %in% c(11)])
}

ymdValues <- ymd(orig$date)
## Remove leap days to keep number of days per season consistent;
## Except for 'season' 5, this just chops off last day of DJF
leapDay <- day(ymdValues) == 29 & month(ymdValues) == 2
orig <- orig[as.numeric(orig$month) %in% season_months[[season]] &
           !leapDay, ]

y <- split(orig$prcp, orig$season_year)
y <- t(sapply(y, function(x) x)) # rows are each year

missing <- is.na(y)
y[missing] <- 0
y <- y/10  # cm
y[y > 20] <- 20 
n <- ncol(y)
nT <- nrow(y)
n_missing <- sum(missing)

maxlag <- 30

if(!allow_missing)
  missing <- matrix(0, nT, n)

if(rounded)  ## inches rounded to nearest 0.01
  y <- round(y/2.54, 2)

if(!rounded && use_inches)
  y <- y/2.54

r <- y

holdout <- matrix(FALSE, nT, n)


doy <- 1:n
dayTime <- paste0(unique(orig$season_year), '-', season_mid[season])
dayIndex <- as.numeric(( ymd(dayTime) - ymd('1850-1-1') ) + 1)
dayCovar <- (dayIndex - mean(dayIndex))/(max(dayIndex) - min(dayIndex))

   

if(FALSE) {
## Starting 1850, through end 2019
    regYr <- 1:365
    leapYr <- c(1:59, 59, 60:365)  # treat Feb 29 as same doy as Feb 28
    nDays <- c(regYr, regYr, rep(c(leapYr, regYr, regYr, regYr), 12), rep(regYr, 4),
               rep(c(leapYr, regYr, regYr, regYr), 29), regYr, regYr, regYr)
    
    dayIndex <- as.numeric(( ymd(orig$time) - ymd('1850-1-1') ) + 1)
    doy <- nDays[dayIndex]
    
    dayCovar <- (dayIndex - mean(dayIndex))/(max(dayIndex) - min(dayIndex))
}
