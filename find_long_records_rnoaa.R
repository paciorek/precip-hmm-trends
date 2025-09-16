## This uses rnoaa package.
## But would still need to clean a la Mark's code.

library(rnoaa)
library(ggplot2)
library(dplyr)

get_cleaned_data <- function(id, date_min = "1899-12-01", date_max = "2021-11-30") {
    out <- meteo_tidy_ghcnd(id, keep_flags = TRUE, var = 'PRCP', date_min = date_min, date_max = date_max)

    if(!nrow(out)) return(NULL)
    out$prcp <- as.numeric(out$prcp)
    
    out$prcp[out$prcp == -9999] <- NA
    ## Values that fail any quality assurance check
    out$prcp[out$qflag_prcp != " "] <- NA
    ## Values with SFLAG = S or missing SFLAG
    out$prcp[out$sflag_prcp == "S" | out$sflag_prcp == " "]  <- NA
    ## Set values with MFLAG = "T" equal to zero
    out$prcp[out$mflag_prcp == "T" ] <- 0
    
    out$prcp <- out$prcp/10 # convert to cm
    return(out)
}
 

if(!exists('stns')) {
    stns <- ghcnd_stations()
    stns <- stns |> dplyr::filter(element == 'PRCP')
    
    save(stns, file = 'ghcnd_stations.Rda')
} else load('ghcnd_stations.Rda')


az_stns <- stns |> dplyr::filter(state == 'AZ')
ca_stns <- stns |> dplyr::filter(state == 'CA')

az_precip <- ghcnd(az_stns$id)
ca_precip <- ghcnd(ca_stns$id)

## rnoaa downloads files to /tmp/.xdg_cache_paciorek/R/noaa_ghcnd

start_year <- 1900
end_year <- 2021
nyr <- end_year - start_year + 1
season_year <- rep(start_year:end_year, each = 31*12)
month <- as.character(rep(rep(c(12,1:11), each = 31), times = nyr))
month[nchar(month) == 1]  <- paste0("0", month[nchar(month) == 1])
day <- as.character(rep(1:31, times = nyr * 12))
day[nchar(day) == 1]  <- paste0("0", day[nchar(day) == 1])
raw_year <- season_year
raw_year[month == 12] <- raw_year[month == 12] - 1 
full_dates <- data.frame(date = paste(raw_year, month, day, sep = "-"), raw_year = raw_year, month = month, season_year = season_year)
exclude <- (month %in% c("04","06","09","11") & day == "31") | (month == "02" & day > "29") |
    (month == "02" & raw_year %% 4 != 0 & day == "29") |
    (month == "02" & raw_year == 1900 & day == "29")
full_dates <- full_dates |> dplyr::filter(!exclude)
full_dates$date <- as.Date(full_dates$date)


az_data <- matrix(0, nrow = nrow(full_dates), ncol = nrow(az_stns))
ca_data <- matrix(0, nrow = nrow(full_dates), ncol = nrow(ca_stns))

## AZ - Phoenix

for(stn in seq_along(az_stns$id)) {
    out <- get_cleaned_data(az_stns$id[stn])
    if(!is.null(out)) {
        out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)
        az_data[ , stn] <- out$prcp
    }
}

save(az_data, az_stns, file = "AZ_ghcnd.Rda")


## PHX

timept <- which(full_dates$date == '1920-1-1')
last <- nrow(full_dates)

dist <- fields::rdist.earth(as.matrix(az_stns[ , c('longitude','latitude')]),
                    matrix(c(-112.1, 33.4), nrow = 1),
                    miles = FALSE)

near <- which(dist < 50)

phx_data <- az_data[ , near]
phx_stns <- az_stns[near, ]

prop <- apply(phx_data[timept:last,], 2, function(x)
    mean(!is.na(x)))
## sub <- ca_precip[ , first < timept & prop > .8]
sub <- phx_data[ , prop > .6] 

time_idx <- seq_len(nrow(phx_data))
first <- apply(sub, 2, function(x) which(!is.na(x))[1])
ord <- order(first)
sub_tmp <- sub[, ord]
plot(time_idx[!is.na(sub_tmp[ , 1])], rep(1, sum(!is.na(sub_tmp[ , 1]))),
     ylim = c(1,ncol(sub_tmp)), xlim = range(time_idx), cex = .4, pch = 16)
for(idx in 2:ncol(sub_tmp))
    points(time_idx[!is.na(sub_tmp[ , idx])], rep(idx, sum(!is.na(sub_tmp[ , idx]))), cex = .4, pch = 16)

## options:

## Litchfield Park good but ends about a year early USC00024977
## Sky Harbor seems good but starts in 1933 USW00023183

## identified Prescott and Winslow as stations analyzed by Zhang with long records and significant p-values

id <- "USC00024977"  # 
out <- get_cleaned_data(id)
out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)

write.csv(out[ , c('date','raw_year','season_year','month','prcp')], file = file.path('data', "AZ_Litchfield.csv"), row.names = FALSE)

## Oroville

dam_locn <- matrix(c(-121.4854488,39.5380359), nrow = 1)
dist <- fields::rdist.earth(as.matrix(ca_stns[ , c('longitude','latitude')]),
                    dam_locn,  miles = FALSE)

near <- which(dist < 125)

oro_stns <- ca_stns[near, ]

plot(oro_stns$longitude, oro_stns$latitude)
points(-121.4854488,39.5380359, col = 'red', pch=16)
points(-121.2695579,40.2988315,col='blue',pch=16)  # Chester
points(-120.404425,40.1641307, col = 'blue',pch=16) # Milford

oro_stns <- oro_stns[oro_stns$latitude > dam_locn[2] & oro_stns$longitude > dam_locn[1],]
oro_data <- matrix(0, nrow = nrow(full_dates), ncol = nrow(oro_stns))

for(stn in seq_along(oro_stns$id)) {
    out <- get_cleaned_data(oro_stns$id[stn])
    out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)
    oro_data[ , stn] <- out$prcp
}

timept <- which(full_dates$date == '1920-1-1')
last <- nrow(full_dates)

prop <- apply(oro_data[timept:last,], 2, function(x)
    mean(!is.na(x)))
## sub <- oro_precip[ , first < timept & prop > .8]
sub <- oro_data[ , prop > .8] 

time_idx <- seq_len(nrow(oro_data))
first <- apply(sub, 2, function(x) which(!is.na(x))[1])
ord <- seq_len(ncol(sub)) # order(first)
sub_tmp <- sub[, ord]
plot(time_idx[!is.na(sub_tmp[ , 1])], rep(1, sum(!is.na(sub_tmp[ , 1]))),
     ylim = c(1,ncol(sub_tmp)), xlim = range(time_idx), cex = .4, pch = 16)
for(idx in 2:ncol(sub_tmp))
    points(time_idx[!is.na(sub_tmp[ , idx])], rep(idx, sum(!is.na(sub_tmp[ , idx]))), cex = .4, pch = 16)


good_stns <- oro_stns[prop > .8,]
plot(good_stns$longitude, good_stns$latitude, xlim = c(-121.5, -120.2), ylim = c(39.5,40.5))
points(-121.4854488,39.5380359, col = 'red', pch=16)
points(-121.2695579,40.2988315,col='blue',pch=16)  # Chester
points(-120.404425,40.1641307, col = 'blue',pch=16) # Milford

## options: about 9 stations with good temporal coverage but most except Quincy are on fringes of watershed


id <- "USC00047195"
out <- get_cleaned_data(id)
out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)

write.csv(out[ , c('date','raw_year','season_year','month','prcp')], file = file.path('data', "CA_Quincy.csv"), row.names = FALSE)


## Bay Area

bay_locn <- matrix(c(-122.3364572,37.8712261), nrow = 1)
dist <- fields::rdist.earth(as.matrix(ca_stns[ , c('longitude','latitude')]),
                    bay_locn,  miles = FALSE)

near <- which(dist < 75)

bay_stns <- ca_stns[near, ]

bay_stns <- bay_stns[bay_stns$latitude > bay_locn[2] & bay_stns$longitude > bay_locn[1],]
bay_data <- matrix(0, nrow = nrow(full_dates), ncol = nrow(bay_stns))

for(stn in seq_along(bay_stns$id)) {
    out <- get_cleaned_data(bay_stns$id[stn])
    out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)
    bay_data[ , stn] <- out$prcp
}

timept <- which(full_dates$date == '1920-1-1')
last <- nrow(full_dates)

prop <- apply(bay_data[timept:last,], 2, function(x)
    mean(!is.na(x)))
## sub <- bay_precip[ , first < timept & prop > .8]
sub <- bay_data[ , prop > .6] 

time_idx <- seq_len(nrow(bay_data))
first <- apply(sub, 2, function(x) which(!is.na(x))[1])
ord <- seq_len(ncol(sub)) # order(first)
sub_tmp <- sub[, ord]
plot(time_idx[!is.na(sub_tmp[ , 1])], rep(1, sum(!is.na(sub_tmp[ , 1]))),
     ylim = c(1,ncol(sub_tmp)), xlim = range(time_idx), cex = .4, pch = 16)
for(idx in 2:ncol(sub_tmp))
    points(time_idx[!is.na(sub_tmp[ , idx])], rep(idx, sum(!is.na(sub_tmp[ , idx]))), cex = .4, pch = 16)

## only Berkeley and Napa State Hosp are options; Vacaville goes back a long time but cuts out early


id <- "USC00040693"
out <- get_cleaned_data(id)
out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)

write.csv(out[ , c('date','raw_year','season_year','month','prcp')], file = file.path('data', "CA_Berkeley.csv"), row.names = FALSE)


## KK Fort Collins, CO station

id <- "USC00053005"
out <- get_cleaned_data(id)
out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)

write.csv(out[ , c('date','raw_year','season_year','month','prcp')], file = file.path('data', "CO_FortCollins.csv"), row.names = FALSE)

## Hazard, KY station

id <- "USC00153714"
out <- get_cleaned_data(id)
out <- merge(out, full_dates, by.x = 'date', by.y = 'date', all.y = TRUE)

write.csv(out[ , c('date','raw_year','season_year','month','prcp')], file = file.path('data', "KY_Hazard.csv"), row.names = FALSE)


## OLD
    
##    dat <- az_precip |> dplyr::filter(id == stn) |> dplyr::filter(element == 'PRCP')
##    value <- dat |> select(starts_with("VALUE"))
##    mflag <- dat |> select(starts_with("MFLAG"))
##    qflag <- dat |> select(starts_with("QFLAG"))
##    sflag <- dat |> select(starts_with("SFLAG"))
        
## Missing values labelled NA
## value[value == -9999] <- NA
## Values that fail any quality assurance check
## value[qflag != " "] <- NA
## Values with SFLAG = S or missing SFLAG
## value[sflag == "S" | sflag == " "]  <- NA
## Set values with MFLAG = "T" equal to zero
## value[mflag == "T" ] <- 0
