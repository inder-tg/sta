
# input: vector, amp.Type, frequency, significance, output: ampCoeffs, slope, linearTrend, pVal
sta.ts <- function(ts, intraAnnualPeriod = c("wetSeason", "drySeason"), adhocPeriod = NULL,
                   ampType = c("mean", "annual", "semi-annual"), freq, 
                   numFreq, delta){
  
  intraAnnualPeriod <- match.arg(intraAnnualPeriod)
  
  ampType <- match.arg(ampType)
  
  if(missing(freq)){
    stop("freq must be provided")
  }

  if(missing(numFreq)){
    stop("numFreq must be provided")
  }

  if(missing(delta)){
    stop("delta must be provided")
  }
  
  if( length(ts) %% freq != 0 ){
    stop("length(ts) must be a multiple of freq")
  }
  
  years <- 1:(length(ts)/freq)
  ampCoeffs <- rep(NA, length(ts)/freq)
  trend <- rep(NA, length(ts)/freq)
  harmFit <- rep(NA, length(ts))
  
  for(i in years){
    fit <- haRmonics(y = ts[( (i-1) * freq + 1 ):( i * freq )], method = "harmR",
                     numFreq = numFreq, delta = delta)
    ampCoeffs[i] <- ifelse(ampType == "mean", fit$amplitude[1], 
                          ifelse(ampType == "annual", fit$amplitude[2], fit$amplitude[3]))
    harmFit[( (i-1) * freq + 1 ):( i * freq )] <- fit$fitted
  }
  
  SENS <- sens.slope(ampCoeffs)
  trend <- median(ampCoeffs - SENS$estimates * 0:(length(years)-1) ) + SENS$estimates * 0:(length(years)-1)
  
  if(is.null(adhocPeriod)){
    basicStats <- getBasicStatsIntra(y = harmFit, freq = freq, intraAnnualPeriod = intraAnnualPeriod)
  } else {
    basicStats <- getBasicStatsAdhoc(y = harmFit, freq = freq, adhocPeriod = adhocPeriod) 
  }
  
  list(ts = ts, fit = harmFit, freq = freq,
       harmCoeffs = ampCoeffs, linearTrend = trend, 
       slope = SENS$estimates, pval = SENS$p.value,
       meanMins = basicStats$min_mean, sdMins = basicStats$min_sd,
       minMins = basicStats$min_min, maxMins = basicStats$min_max,
       meanMaxs = basicStats$max_mean, sdMaxs = basicStats$max_sd,
       minMaxs = basicStats$max_min, maxMaxs = basicStats$max_max)
}

# in order to speed things up we are gonna run our code (sta.ts) in parallel
# input: df, amp.Type, frequency, significance, parallel, noCores, save
# output: ampCoeffs, slope, linearTrend, pVal
globalVariables('i')
sta.matrix <- function(df,
                       intraAnnualPeriod = c("wetSeason", "drySeason"), 
                       adhocPeriod = NULL, freq, numFreq, delta, 
                       numCores = 20, dirToSaveSTA = NULL){
  
  if( !is.null(adhocPeriod) ){
    intraAnnualPeriod <- NULL
  }
  
  NROW <- nrow(df)
  NCOL <- ncol(df)
  mean <- matrix(NA, nrow = nrow(df), ncol = 4) # X,Y,estimate,pVal
  mean[,c(1:2)] <- df[, c(1:2)]
  
  mean_basicStats <- matrix(NA, nrow = nrow(df), ncol = 10) # X,Y,estimate,pVal
  mean_basicStats[,c(1:2)] <- df[, c(1:2)]
  
  annual <- matrix(NA, nrow = nrow(df), ncol = 4) # X,Y,estimate,pVal
  annual[,c(1:2)] <- df[, c(1:2)]
  
  annual_basicStats <- matrix(NA, nrow = nrow(df), ncol = 10) # X,Y,estimate,pVal
  annual_basicStats[,c(1:2)] <- df[, c(1:2)]
  
  semiannual <- matrix(NA, nrow = nrow(df), ncol = 4) # X,Y,estimate,pVal
  semiannual[,c(1:2)] <- df[, c(1:2)]
  
  semiannual_basicStats <- matrix(NA, nrow = nrow(df), ncol = 10) # X,Y,estimate,pVal
  semiannual_basicStats[,c(1:2)] <- df[, c(1:2)]
  
  if( !is.null(dirToSaveSTA) ){
    reportFileName <- paste0( dirToSaveSTA, "/sta_progress.txt" )
    file.create(path = reportFileName, showWarnings = FALSE)
    write( "--- STA began at ---", file = reportFileName, append = TRUE )
    write( as.character(Sys.time()[1]), file = reportFileName, append = TRUE )
  }
  
  kluster <- parallel::makeCluster(spec = numCores, outfile = "")
  registerDoParallel(kluster)
  
  output <- foreach(i = 1:NROW, .combine = "rbind", .export = c("sta.ts", "getBasicStatsAdhoc", "getBasicStatsIntra", "getDaysToAnalyze"),
                    .packages = c("raster", "geoTS", "trend") ) %dopar% {
                      getSTA_mean <- sta.ts(df[i,3:NCOL],
                                            intraAnnualPeriod = intraAnnualPeriod, 
                                            adhocPeriod = adhocPeriod,
                                            ampType = "mean", freq = freq, 
                                            numFreq = numFreq, delta = delta)
                      
                       getSTA_annual <- sta.ts(df[i,3:NCOL],
                                              intraAnnualPeriod = intraAnnualPeriod, 
                                              adhocPeriod = adhocPeriod,
                                              ampType = "annual", freq = freq,
                                              numFreq = numFreq, delta = delta)
                      
                      getSTA_semiannual <- sta.ts(df[i,3:NCOL],
                                                  intraAnnualPeriod = intraAnnualPeriod, 
                                                  adhocPeriod = adhocPeriod,
                                                  ampType = "semi-annual", freq = freq, 
                                                  numFreq = numFreq, delta = delta)
                      
                      s <- c(getSTA_mean$slope, getSTA_mean$pval,
                             getSTA_mean$meanMins, getSTA_mean$sdMins,
                             getSTA_mean$minMins, getSTA_mean$maxMins,
                             getSTA_mean$meanMaxs, getSTA_mean$sdMaxs,
                             getSTA_mean$minMaxs, getSTA_mean$maxMaxs,
                             getSTA_annual$slope, getSTA_annual$pval,
                             getSTA_annual$meanMins, getSTA_annual$sdMins,
                             getSTA_annual$minMins, getSTA_annual$maxMins,
                             getSTA_annual$meanMaxs, getSTA_annual$sdMaxs,
                             getSTA_annual$minMaxs, getSTA_annual$maxMaxs,
                             getSTA_semiannual$slope, getSTA_semiannual$pval,
                             getSTA_semiannual$meanMins, getSTA_semiannual$sdMins,
                             getSTA_semiannual$minMins, getSTA_semiannual$maxMins,
                             getSTA_semiannual$meanMaxs, getSTA_semiannual$sdMaxs,
                             getSTA_semiannual$minMaxs, getSTA_semiannual$maxMaxs)
                      
                      if( !is.null(dirToSaveSTA) ){
                        texto <- paste0("Working on ROW:", i)
                        write(texto, file = reportFileName, append = T)
                      }
                      
                      return(s)
                    }
  stopCluster(kluster)
  if( !is.null(dirToSaveSTA) ){
    write( as.character(Sys.time()[1]), file = reportFileName, append = TRUE)
    write("If save = TRUE, output will be saved in this directory shortly", file = reportFileName, append = TRUE)
    write( "--- STA ended at ---", file = reportFileName, append = TRUE)
  }
  
  mean[, c(3,4)] <- output[,1:2]
  mean_basicStats[, 3:10] <- output[,3:10]
  #---
  annual[, c(3,4)] <- output[,11:12]
  annual_basicStats[, 3:10] <- output[,13:20]
  #---
  semiannual[, c(3,4)] <- output[,21:22]
  semiannual_basicStats[,3:10] <- output[,23:30]
  
  list(mean = mean, mean_basicStats = mean_basicStats,
       annual = annual, annual_basicStats = annual_basicStats,
       semiannual = semiannual, semiannual_basicStats = semiannual_basicStats)
}

# input: data, interAnnualPeriod, output: data chopped
getData <- function(data, freq, startYear, endYear, interAnnualPeriod){
  output <- NULL
  
  beginPeriod <- (interAnnualPeriod-1) * freq
  endPeriod <- interAnnualPeriod * freq

  yearsToAnalyze <- c(sapply(1:length(beginPeriod), function(s) (beginPeriod[s]+1):endPeriod[s] ))
  years <- (startYear:endYear)[interAnnualPeriod]
  
  if( inherits(data, "numeric") ){# class(data) == "numeric"
    output <- data[yearsToAnalyze]
  }
  
  if( inherits(data, "matrix") ){ # class(data) == "matrix"
    output <- data[,yearsToAnalyze+2]
  }
  
  list(data = output, analyzedPeriod = yearsToAnalyze, years = years)
}

getDaysToAnalyze <- function(intraAnnualPeriod, years, freq){
  
  if( !(freq == 12 | freq == 23 | freq == 36) ){
    stop("freq must be 12, 23 or 36")
  }
  
  daysToAnalyze <- list()
  
  if(freq == 36){
    if(intraAnnualPeriod == "wetSeason"){
      wetSeasonBegin <- 13 + 0:(length(years)-1) * freq
      wetSeasonEnd <- 30 + 0:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(wetSeasonBegin), function(s) c(wetSeasonBegin[s],wetSeasonEnd[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(wetSeasonBegin), function(s) wetSeasonBegin[s]:wetSeasonEnd[s] ))
    }
    
    if(intraAnnualPeriod == "drySeason"){
      drySeasonBeginA <- 31 + 0:(length(years)-2) * freq
      drySeasonEndA <- 36 + 0:(length(years)-2) * freq
      drySeasonBeginB <- 1 + 1:(length(years)-1) * freq
      drySeasonEndB <- 12 + 1:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(drySeasonBeginA), 
                                        function(s) c(drySeasonBeginA[s], drySeasonEndB[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(drySeasonBeginA),
                                function(s) c(drySeasonBeginA[s]:drySeasonEndA[s],
                                              drySeasonBeginB[s]:drySeasonEndB[s]) ))
    }
  }
  
  if(freq == 12){
    if(intraAnnualPeriod == "wetSeason"){
      wetSeasonBegin <- 5 + 0:(length(years)-1) * freq
      wetSeasonEnd <- 10 + 0:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(wetSeasonBegin), function(s) c(wetSeasonBegin[s], wetSeasonEnd[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(wetSeasonBegin), function(s) wetSeasonBegin[s]:wetSeasonEnd[s] ))
    }
    
    if(intraAnnualPeriod == "drySeason"){
      drySeasonBeginA <- 11 + 0:(length(years)-2) * freq
      drySeasonEndA <- 12 + 0:(length(years)-2) * freq
      drySeasonBeginB <- 1 + 1:(length(years)-1) * freq
      drySeasonEndB <- 4 + 1:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(drySeasonBeginA), 
                                        function(s) c(drySeasonBeginA[s], drySeasonEndB[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(drySeasonBeginA),
                                function(s) c(drySeasonBeginA[s]:drySeasonEndA[s],
                                              drySeasonBeginB[s]:drySeasonEndB[s]) ))
    }
  }
  
  if(freq == 23){
    if(intraAnnualPeriod == "wetSeason"){
      wetSeasonBegin <- 9 + 0:(length(years)-1) * freq
      wetSeasonEnd <- 19 + 0:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(wetSeasonBegin), function(s) c(wetSeasonBegin[s], wetSeasonEnd[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(wetSeasonBegin), function(s) wetSeasonBegin[s]:wetSeasonEnd[s] ))
    }
    
    if(intraAnnualPeriod == "drySeason"){
      drySeasonBeginA <- 20 + 0:(length(years)-2) * freq
      drySeasonEndA <- 23 + 0:(length(years)-2) * freq
      drySeasonBeginB <- 1 + 1:(length(years)-1) * freq
      drySeasonEndB <- 8 + 1:(length(years)-1) * freq
      daysToAnalyze$partial <- c(sapply(1:length(drySeasonBeginA), 
                                function(s) c(drySeasonBeginA[s], drySeasonEndB[s]) ))
      daysToAnalyze$full <- c(sapply(1:length(drySeasonBeginA),
                                     function(s) c(drySeasonBeginA[s]:drySeasonEndA[s],
                                                   drySeasonBeginB[s]:drySeasonEndB[s]) ))
    }
  }
  
  daysToAnalyze
}

getGlobalDaysUsedBasicStats <- function(intraAnnualPeriod = c("wetSeason", "drySeason"), 
                                        adhocPeriod = NULL, 
                                        globalDays, years, freq){
  
  intraAnnualPeriod <- match.arg(intraAnnualPeriod)
  
  if(!is.null(adhocPeriod)){
    
    if( !inherits(adhocPeriod, "list") | length(adhocPeriod) != 2 ){ # class(adhocPeriod) != "list"
      stop("adhocPeriod must be a list of length 2")
    } else {
      daysToAnalyze <- adhocPeriod$partial
    }
    
  } else {
    daysToAnalyze <- getDaysToAnalyze(intraAnnualPeriod = intraAnnualPeriod, 
                                      years = years, freq = freq)$partial
  }
  
  globalDays[daysToAnalyze]
}

getBasicStatsIntra <- function(y, freq, intraAnnualPeriod){
  years <- 1:(length(y)/freq)
  daysToAnalyze <- getDaysToAnalyze(intraAnnualPeriod = intraAnnualPeriod, 
                                    years = years, freq = freq)$partial
  
  if(intraAnnualPeriod == "wetSeason"){
    mins <- sapply(seq(1,2*length(years),2), function(s) min(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
    maxs <- sapply(seq(1,2*length(years),2), function(s) max(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
    mean_mins <- mean(mins, na.rm = T)
    sd_mins <- sd(mins, na.rm = T)
    min_mins <- min(mins, na.rm = T)
    max_mins <- max(mins, na.rm = T)
    mean_maxs <- mean(maxs, na.rm = T)
    sd_maxs <- sd(maxs, na.rm = T)
    min_maxs <- min(maxs, na.rm = T)
    max_maxs <- max(maxs, na.rm = T)
  }
  
  if(intraAnnualPeriod == "drySeason"){
    mins <- sapply(seq(1,2*length(years)-2,2), function(s) min(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
    maxs <- sapply(seq(1,2*length(years)-2,2), function(s) max(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
    mean_mins <- mean(mins, na.rm = T)
    sd_mins <- sd(mins, na.rm = T)
    min_mins <- min(mins, na.rm = T)
    max_mins <- max(mins, na.rm = T)
    mean_maxs <- mean(maxs, na.rm = T)
    sd_maxs <- sd(maxs, na.rm = T)
    min_maxs <- min(maxs, na.rm = T)
    max_maxs <- max(maxs, na.rm = T)
  }
  
  list(min_mean = mean_mins, min_sd = sd_mins, min_min = min_mins, min_max = max_mins,
       max_mean = mean_maxs, max_sd = sd_maxs, max_min = min_maxs, max_max = max_maxs,
       daysToAnalyze = daysToAnalyze)
}

getBasicStatsAdhoc <- function(y, freq, adhocPeriod){
  
  years <- 1:(length(y)/freq)
  
  if( !inherits(adhocPeriod, "list") | length(adhocPeriod) != 2 ){ # class(adhocPeriod) != "list"
    stop("adhocPeriod must be a list of length 2")
  }
  
  daysToAnalyze <- adhocPeriod$partial
  mins <- sapply(seq(1,length(daysToAnalyze),2), function(s) min(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
  maxs <- sapply(seq(1,length(daysToAnalyze),2), function(s) max(y[daysToAnalyze[s]:daysToAnalyze[s+1]] ))
  mean_mins <- mean(mins, na.rm = T)
  sd_mins <- sd(mins, na.rm = T)
  min_mins <- min(mins, na.rm = T)
  max_mins <- max(mins, na.rm = T)
  mean_maxs <- mean(maxs, na.rm = T)
  sd_maxs <- sd(maxs, na.rm = T)
  min_maxs <- min(maxs, na.rm = T)
  max_maxs <- max(maxs, na.rm = T)
  
  list(min_mean = mean_mins, min_sd = sd_mins, min_min = min_mins, min_max = max_mins,
       max_mean = mean_maxs, max_sd = sd_maxs, max_min = min_maxs, max_max = max_maxs,
       daysToAnalyze = daysToAnalyze)
}

getOutput <- function(data,
                      intraAnnualPeriod, adhocPeriod,
                      freq, numFreq, delta, numCores, dirToSaveSTA){
  output <- list()
  
  # -- changing to intraAnnualPeriod, adhocPeriod = NULL
  
  if(!is.null(adhocPeriod)){
    intraAnnualPeriod <- NULL
  } 
  
  if( inherits(data, "numeric") ){ # class(data) == "numeric"
    output$mean <- sta.ts(ts = data,
                          intraAnnualPeriod = intraAnnualPeriod, adhocPeriod = adhocPeriod,
                          ampType = "mean", freq = freq, numFreq = numFreq,
                          delta = delta)

    output$annual <- sta.ts(ts = data,
                            intraAnnualPeriod = intraAnnualPeriod, adhocPeriod = adhocPeriod,
                            ampType = "annual", freq = freq, numFreq = numFreq, delta = delta)

    output$semiannual <- sta.ts(ts = data,
                                intraAnnualPeriod = intraAnnualPeriod, adhocPeriod = adhocPeriod,
                                ampType = "semi-annual", freq = freq, numFreq = numFreq, delta = delta)
  }
  
  if( inherits(data, "matrix")){ #  class(data) == "matrix"
    output <- sta.matrix(df = data,
                         intraAnnualPeriod = intraAnnualPeriod, adhocPeriod = adhocPeriod,
                         freq = freq, numFreq = numFreq, delta = delta, numCores = numCores,
                         dirToSaveSTA = dirToSaveSTA)
  }
  
  output
}

get.sta.numeric <- function(data, days,
                            freq, numFreq = 4, delta = 0.2, 
                            startYear = 2000, endYear = 2018,
                            intraAnnualPeriod = c("wetSeason", "drySeason"), 
                            interAnnualPeriod, adhocPeriod = NULL){
  
  output <- list()
  output$data <- data # backUp
  output$freq <- freq
  output$startYear <- startYear
  output$endYear <- endYear
  
  intraAnnualPeriod <- match.arg(intraAnnualPeriod)
  output$intraAnnualPeriod <- intraAnnualPeriod
  output$interAnnualPeriod <- interAnnualPeriod
  
  data <- getData(data = data, freq = freq, startYear = startYear, endYear = endYear, 
                  interAnnualPeriod = interAnnualPeriod)
  
  # globalDaysToAnalyze <- getData(data = as.numeric(days), freq = freq, interAnnualPeriod = interAnnualPeriod)
  globalDaysToAnalyze <- getData(data = days, freq = freq, startYear = startYear, endYear = endYear, 
                                 interAnnualPeriod = interAnnualPeriod)
  
  aux_output <- list()
  
  if(is.null(adhocPeriod)){
    aux_output$intervalsUsedBasicStats <- getGlobalDaysUsedBasicStats(intraAnnualPeriod = intraAnnualPeriod,
                                                                      globalDays = globalDaysToAnalyze$data,
                                                                      adhocPeriod = NULL, years = data$years,
                                                                      freq = freq)
    
    aux_output$sta <- getOutput(data = data$data,
                                intraAnnualPeriod = intraAnnualPeriod, adhocPeriod = NULL,
                                freq = freq, numFreq = numFreq, delta = delta)
  } else {
    aux_output$intervalsUsedBasicStats <- getGlobalDaysUsedBasicStats(adhocPeriod = adhocPeriod,
                                                                      # globalDays = globalDaysToAnalyze$data,
                                                                      globalDays = 1:(length(startYear:endYear) * freq),
                                                                      years = data$years, freq = freq)
    
    aux_output$sta <- getOutput(data = data$data,
                                adhocPeriod = adhocPeriod,
                                freq = freq, numFreq = numFreq, delta = delta)
  }
  
  output$daysUsedFit <- data$analyzedPeriod
  output$intervalsUsedBasicStats <- aux_output$intervalsUsedBasicStats
  
  output$ts <- aux_output$sta$mean$ts
  output$fit <- aux_output$sta$mean$fit
  
  aux_output$sta$mean$freq <- NULL
  aux_output$sta$mean$ts <- NULL
  aux_output$sta$mean$fit <- NULL
  
  aux_output$sta$annual$freq <- NULL
  aux_output$sta$annual$ts <- NULL
  aux_output$sta$annual$fit <- NULL
  
  aux_output$sta$semiannual$freq <- NULL
  aux_output$sta$semiannual$ts <- NULL
  aux_output$sta$semiannual$fit <- NULL
  
  output$sta <- aux_output$sta
  
  output
}

get.sta.matrix <- function(data, days, freq, numFreq = 4, delta = 0.2,
                           startYear = startYear, endYear = endYear,
                           intraAnnualPeriod = c("wetSeason", "drySeason"),
                           interAnnualPeriod, adhocPeriod = NULL, 
                           numCores = 20, dirToSaveSTA = NULL){
  
  dataAux <- getData(data = data, freq = freq, startYear = startYear, endYear = endYear,
                     interAnnualPeriod = interAnnualPeriod)
  
  globalDaysToAnalyze <- getData(data = days, freq = freq, startYear = startYear, endYear = endYear,
                                 interAnnualPeriod = interAnnualPeriod)
  
  aux_output <- list()
  if(is.null(adhocPeriod)){
    aux_output$intervalsUsedBasicStats <- getGlobalDaysUsedBasicStats(intraAnnualPeriod = intraAnnualPeriod,
                                                                      globalDays = globalDaysToAnalyze$data,
                                                                      adhocPeriod = NULL, years = dataAux$years,
                                                                      freq = freq)
    aux_output$sta <- getOutput(data = data, intraAnnualPeriod = intraAnnualPeriod,
                                adhocPeriod = NULL, freq = freq, numFreq = numFreq, delta = delta,
                                dirToSaveSTA = dirToSaveSTA, numCores = numCores)
  } else {
    aux_output$intervalsUsedBasicStats <- getGlobalDaysUsedBasicStats(adhocPeriod = adhocPeriod,
                                                                      # globalDays = globalDaysToAnalyze$data,
                                                                      globalDays = 1:(length(startYear:endYear) * freq),
                                                                      years = dataAux$years, freq = freq)
    aux_output$sta <- getOutput(data = data, adhocPeriod = adhocPeriod,
                                freq = freq, numFreq = numFreq, delta = delta, 
                                dirToSaveSTA = dirToSaveSTA, numCores = numCores)
  }
  
  output <- list()
  output$freq <- freq
  output$daysUsedFit <- dataAux$analyzedPeriod # previouls data$analyzedPeriod
  output$intervalsUsedBasicStats <- aux_output$intervalsUsedBasicStats
  output$sta <- aux_output$sta
  
  output
}

getMapview <- function(mapRaster, colPal, nameLayer, 
                       typeQuery = c("mousemove", "click"), 
                       label = T){ 

  palTest <- colorRampPalette( colors = colPal, alpha = T )
  palTest_rev <- colorRampPalette( colors = rev(colPal) )
  
  typeQuery <- match.arg(typeQuery)
  
  m <- mapview(mapRaster, na.color = "transparent", 
               col.regions = palTest(256), 
               map.types = c("Esri.WorldImagery", 
                             "CartoDB.Positron", 
                             "OpenStreetMap", 
                             "OpenTopoMap"), 
               query.type = typeQuery,
               label = label, homebutton = F,
               layer.name = nameLayer, legend = T)
  
  m_rev <- mapview(mapRaster, na.color = "transparent", 
                   col.regions = palTest_rev(256), 
                   map.types = c("Esri.WorldImagery", 
                                 "OpenTopoMap"), 
                   homebutton = F,
                   layer.name = nameLayer, legend = T)
  
  m@map[[1]]$calls[[8]]$args[[1]]$labels <- rev(m@map[[1]]$calls[[8]]$args[[1]]$labels)
  
  m@map[[1]]$calls[[8]]$args[[1]]$title = paste0("slope-", nameLayer) #"Slope"
  
  m@map[[1]]$calls[[8]]$args[[1]]$colors <- m_rev@map[[1]]$calls[[6]]$args[[1]]$colors
  
  m
}

matrixToRaster <- function(matrix, raster){
  matrix_TEMP <- data.frame(x = matrix[,1], y = matrix[,2], values = matrix[,3])
  
  sp::coordinates(matrix_TEMP) <- ~x + y
  
  sp::gridded(matrix_TEMP) <- TRUE
  
  raster_TEMP <- raster(matrix_TEMP)
  
  raster::projection(raster_TEMP) <- raster::projection(raster)
  
  raster_TEMP
}

matrixToMapview <- function(mat, significance, master, 
                            color = c("Reds", "Greens", "Blues"), 
                            label){
  
  color <- match.arg(color)
  
  if( !is.null(significance) ){
    mat[,3][ mat[,4] >= significance ] <- NA
  }

  RASTER <- matrixToRaster(matrix = mat, raster = master)
  getMapview(mapRaster = RASTER, colPal = brewer.pal(256, color), 
             nameLayer = label, typeQuery = "click")
  
}

saveMaps <- function(output, master, dirToSaveSTA){
  
  raster_mean_slope <- matrixToRaster(matrix = output$sta$mean[,c(1,2,3)], 
                                      raster = master)
  raster_mean_pval <- matrixToRaster(matrix = output$sta$mean[,c(1,2,4)], 
                                     raster = master)
  raster_mean_meanMins <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,3)], 
                                         raster = master)
  raster_mean_sdMins <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,4)], 
                                       raster = master)
  raster_mean_minMins <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,5)], 
                                        raster = master)
  raster_mean_maxMins <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,6)], 
                                        raster = master)
  raster_mean_meanMaxs <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,7)], 
                                         raster = master)
  raster_mean_sdMaxs <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,8)], 
                                       raster = master)
  raster_mean_minMaxs <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,9)], 
                                        raster = master)
  raster_mean_maxMaxs <- matrixToRaster(matrix = output$sta$mean_basicStats[,c(1,2,10)], 
                                        raster = master)
  writeRaster(raster_mean_slope, 
              filename = paste0(dirToSaveSTA, "/sta_mean_slope"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_pval, 
              filename = paste0(dirToSaveSTA, "/sta_mean_pval"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_meanMins, 
              filename = paste0(dirToSaveSTA, "/sta_mean_meanMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_sdMins, 
              filename = paste0(dirToSaveSTA, "/sta_mean_sdMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_minMins, 
              filename = paste0(dirToSaveSTA, "/sta_mean_minMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_maxMins, 
              filename = paste0(dirToSaveSTA, "/sta_mean_maxMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_meanMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_mean_meanMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_sdMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_mean_sdMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_minMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_mean_minMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_mean_maxMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_mean_maxMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  
  # ---
  
  raster_annual_slope <- matrixToRaster(matrix = output$sta$annual[,c(1,2,3)], 
                                      raster = master)
  raster_annual_pval <- matrixToRaster(matrix = output$sta$annual[,c(1,2,4)], 
                                     raster = master)
  raster_annual_meanMins <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,3)], 
                                         raster = master)
  raster_annual_sdMins <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,4)], 
                                       raster = master)
  raster_annual_minMins <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,5)], 
                                        raster = master)
  raster_annual_maxMins <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,6)], 
                                        raster = master)
  raster_annual_meanMaxs <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,7)], 
                                         raster = master)
  raster_annual_sdMaxs <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,8)], 
                                       raster = master)
  raster_annual_minMaxs <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,9)], 
                                        raster = master)
  raster_annual_maxMaxs <- matrixToRaster(matrix = output$sta$annual_basicStats[,c(1,2,10)], 
                                        raster = master)
  writeRaster(raster_annual_slope, 
              filename = paste0(dirToSaveSTA, "/sta_annual_slope"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_pval, 
              filename = paste0(dirToSaveSTA, "/sta_annual_pval"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_meanMins, 
              filename = paste0(dirToSaveSTA, "/sta_annual_meanMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_sdMins, 
              filename = paste0(dirToSaveSTA, "/sta_annual_sdMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_minMins, 
              filename = paste0(dirToSaveSTA, "/sta_annual_minMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_maxMins, 
              filename = paste0(dirToSaveSTA, "/sta_annual_maxMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_meanMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_annual_meanMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_sdMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_annual_sdMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_minMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_annual_minMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_annual_maxMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_annual_maxMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  
  # ---
  
  raster_semiannual_slope <- matrixToRaster(matrix = output$sta$semiannual[,c(1,2,3)], 
                                        raster = master)
  raster_semiannual_pval <- matrixToRaster(matrix = output$sta$semiannual[,c(1,2,4)], 
                                       raster = master)
  raster_semiannual_meanMins <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,3)], 
                                           raster = master)
  raster_semiannual_sdMins <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,4)], 
                                         raster = master)
  raster_semiannual_minMins <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,5)], 
                                          raster = master)
  raster_semiannual_maxMins <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,6)], 
                                          raster = master)
  raster_semiannual_meanMaxs <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,7)], 
                                           raster = master)
  raster_semiannual_sdMaxs <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,8)], 
                                         raster = master)
  raster_semiannual_minMaxs <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,9)], 
                                          raster = master)
  raster_semiannual_maxMaxs <- matrixToRaster(matrix = output$sta$semiannual_basicStats[,c(1,2,10)], 
                                          raster = master)
  writeRaster(raster_semiannual_slope, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_slope"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_pval, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_pval"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_meanMins, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_meanMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_sdMins, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_sdMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_minMins, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_minMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_maxMins, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_maxMins"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_meanMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_meanMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_sdMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_sdMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_minMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_minMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
  writeRaster(raster_semiannual_maxMaxs, 
              filename = paste0(dirToSaveSTA, "/sta_semiannual_maxMaxs"),
              format = "GTiff", datatype = "FLT4S",
              overwrite = T)
}

#' Get a \code{RasterLayer} with no missing values from a \code{Raster*} object
#' 
#' The term \emph{master} refers to a raster layer whose extent and coordinate reference system
#' are used as a reference to rasterize further objects, e.g. matrices. To rasterize, \emph{master} must be free of missing values.
#' 
#' @param x \code{Raster*} object
#' 
#' @export
#' 
#' @importFrom raster subset
#' @importFrom raster nlayers
#' 
#' @seealso \code{\link[geoTS]{matrixToRaster}}
#'           
#' @return \code{RasterLayer}
getMaster <- function(x){
  LAYER <- subset(x, sample(1:nlayers(x),1))
  LAYER[is.na(LAYER)] <- 1e4
  LAYER
}
