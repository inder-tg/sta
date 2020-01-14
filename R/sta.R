#' Statistical Seasonal Trend Analysis
#' 
#' @param data              numeric vector, matrix or RasterStack object
#' @param freq              integer with the number of observations per period. See \code{Details}
#' @param numFreq           integer with the number of frequencies to employ in harmonic regression fitting. 
#'                          See \code{\link[geoTS]{haRmonics}}
#' @param delta             numeric (positive) controlling regularization and prevent non-invertible
#'                          hat matrix in harmonic regression model. See \code{\link[geoTS]{haRmonics}}
#' @param startYear         numeric, time series initial year                            
#' @param endYear           numeric, time series last year                            
#' @param intraAnnualPeriod character indicating seasons (wet or dry) to be considered for additional 
#'                          statistical analysis. See \code{Details}
#' @param interAnnualPeriod numeric vector indicating the number of years to be considered in STA. For instance,
#'                          1:5, indicates that the first five years will be
#'                          utilized for STA. Similarly, c(2,6,10) indicates that the second, sixth and tenth
#'                          years will be utilized for STA. See \code{Details}
#' @param adhocPeriod       numeric vector with the specific observations to be considered in additional
#'                          statistical analysis. See \code{Details}
#' @param plot              logical, should STA maps be plotted, default is \code{FALSE}. See \code{Details}
#' @param significance      numeric, significance level used in additional statistical analysis
#' @param save              logical, should STA output be saved, default is \code{FALSE}
#' @param dirToSaveSTA      character with full path name, where the STA maps will be saved, applicable when 
#'                          \code{save = TRUE}
#' @param numCores          integer given the number of cores to use; pertinent when \code{data} is a
#'                          \code{RasterStack} or a \code{matrix}
#' @param master            character with full path name of a \emph{.tif} file whose extend and \code{\link[raster]{crs}}
#'                          are used to rasterize \code{sta} output. See \code{Details}
#' 
#' @export
#' 
#' @examples
#' # Use of interAnnualPeriod
#' sta_wetSeason_21016 <- sta(data = wetlands, freq = 36, interAnnualPeriod = c(2, 10, 16), 
#'                            plot = TRUE)
#' str(sta_wetSeason_21016)
#' 
#' # Use of adhocPeriod; include all the observations for calculating basic statistics
#' adhoc <- list()
#' beginPeriod <- (1:19-1) * 36
#' endPeriod <- 1:19 * 36 # adhoc$partial length must be even
#' adhoc$partial <- c( sapply(1:length(beginPeriod), function(s) c(beginPeriod[s]+1, endPeriod[s]) ) )
#' adhoc$full <- c( sapply(1:length(beginPeriod), function(s) (beginPeriod[s]+1):endPeriod[s]) )
#' sta(data = wetlands, freq = 36, adhocPeriod = adhoc, plot = TRUE, significance = 0.05)
#' 
#' @importFrom raster rasterToPoints
#' @importFrom raster writeRaster
#' @importFrom raster subset
#' @importFrom raster stack
#' @importFrom raster nlayers
#' @importFrom graphics layout
#' @importFrom graphics plot
#' @importFrom graphics axis
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics title
#' @importFrom stats ts
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom geoTS haRmonics
#' @importFrom trend sens.slope
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom utils globalVariables
#' @importFrom mapview mapview
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' 
#' @details When the input is a \code{\link[base]{matrix}}, its first two columns must correspond 
#' to geographic coordinates. For instance, the matrix resulting from applying \code{\link[raster]{rasterToPoints}} 
#' to a \emph{RasterStack} follows this rule. 
#' 
#' \code{freq} must be either 12 (monthly observations), 23 (Landsat annual scale) or 36 (10-day composite)
#' as this version implements STA for time series with these frequencies.
#' 
#' This version sets \code{intraAnnualPeriod} to either the \code{wetSeason} or the \code{drySeason} 
#' of Mexico. Empirical evidence suggests that while wet season runs from May to October, dry season
#' runs from November to April. Should a desired STA require specific months/days, these must be
#' provided through \code{adhocPeriod}.
#' 
#' When \code{interAnnualPeriod} is not specified and \code{class(data)=numeric},
#' \code{interAnnualPeriod = 1:(length(data)/freq)}; when \code{class(data)} is either \code{RasterStack} 
#' or \code{matrix}, \code{interAnnualPeriod = 1:((ncol(data)-2)/freq)}.
#' 
#' Since \code{adhocPeriod} defines an inter annual period "ad-hoc", the specific days of this ad-hoc
#' season must be known in advance and consequently, the specific time-points (with respect to the 
#' time series under consideration) must be provided in a numeric vector.
#' 
#' When \code{plot=T} and \code{class(data)=numeric}, a generic plot is displayed. The plot area is
#' divided in four rectangles of equal dimensions using \code{\link[graphics]{layout}}. The rectangle in
#' the top displays the original time series (gray points) along with the harmonic regression fit (black line)
#' based on the observations indicated by \code{interAnnualPeriod}, also there are highlighted those
#' points (red) given by \code{intraAnnualPeriod} or \code{adhocPeriod} which are used to calculate
#' basic statistics. 
#' In the rectangle below
#' the output of the statistical trend analysis (performed with \code{\link[trend]{sens.slope}}) on the 
#' shape parameter mean is displayed; for this analysis only the years given by \code{interAnnualPeriod} 
#' are considered; a p-value smaller than \code{significance} is highlighted by displaying this value in red.
#' In the following two rectangles, similar STA (for annual and semi-annual shape parameters) are displayed.
#' 
#' When \code{plot=T} and \code{class(data)} is RasterStack or matrix, a \code{\link[mapview]{mapview}}
#' object is plotted. This mapview object contains the slope analysis maps of the mean, annual and semi-annual
#' shape parameters. 
#' 
#' When \code{save=T}, a valid \code{dirToSaveSTA} must be provided, that is, this folder should have been
#' created previously. In this case, the STA output is saved in \code{dirToSaveSTA}. This version
#' saves arrays of STA of the mean, annual and semi-annual parameters (along with their corresponding basic statistics)
#' in the file \code{sta_matrix_output.RData} inside \code{dirToSaveSTA}. Also, in the same directory,
#' the file \code{sta_progress.txt} records the progress of the STA process.
#' 
#' \code{master} provides the filename (with complete path) of an image free of missing values (\code{NA}).
#' The extent of this image must coincide with that of the RasterStack under analysis.
#' 
#' @note \code{save=T}, \code{dirToSaveSTA}, \code{numCores} and \code{master} are required when \code{data} is either a
#' RasterStack or a matrix. The aforementioned \emph{basic statistics} are: mean and standard deviation
#' of the time series of annual maximum and minimum as well as the global minima and maxima.
#' 
#' @return A list containing a bunch of elements
#' 
#' @references Eastman, R., Sangermano, F., Ghimine, B., Zhu, H., Chen, H., Neeti, N., Cai, Y., Machado E., Crema, S. (2009).
#' \emph{Seasonal trend analysis of image time series},
#' International Journal of Remote Sensing, \bold{30(10)}, 2721--2726.
#' 
# sta <- function(data, 
#                 freq, numFreq = 4, delta = 0.2, 
#                 startYear = 2000, endYear = 2019,
#                 intraAnnualPeriod = c("wetSeason", "drySeason"), 
#                 interAnnualPeriod, adhocPeriod = NULL,
#                 plot = F, significance = NULL){
sta <- function(data, 
                freq, numFreq = 4, delta = 0.2, 
                startYear = 2000, endYear = 2019,
                intraAnnualPeriod = c("wetSeason", "drySeason"), 
                interAnnualPeriod, adhocPeriod = NULL,
                plot = F, significance = NULL,
                save = F, dirToSaveSTA = NULL, numCores = 20, master = NULL){

  if( !(class(data) == "numeric" | class(data) == "RasterStack" | class(data) == "matrix") ){
    stop("'data' must be either numeric, RasterStack or matrix")
  }
  
  if( class(data) == "RasterStack" ){
    if( missing(master) ){
      message("'master' was not provided, calculating...")
      master <- getMaster(data)
    }
    message("Transferring 'data' to matrix...")
    data <- rasterToPoints(data)
  }
  
  if( class(data) == "matrix" & missing( master ) ){
    stop("'data' is of matrix class, hence 'master' must be provided")
  }
  
  if( missing(freq) ){
    stop("'freq' must be provided")
  }
  
  if( missing(interAnnualPeriod) ){
    message("Calculating 'interAnnualPeriod'...")
    if( class(data) == "numeric" ){
      interAnnualPeriod <- 1:(length(data)/freq)
    } else {
      interAnnualPeriod <- 1:((ncol(data)-2)/freq)
    }
  }
  
  if( !is.null(adhocPeriod) ){
    if( !is.list(adhocPeriod) | length(adhocPeriod) != 2 ){
      stop("'adhocPeriod' must be a list of lenght 2")
    }
  }
  
  if( class(data) == "numeric" ){
    message("Calculating output for class numeric...")
    days <- as.numeric(1:length(data))
    output <- get.sta.numeric(data = data, days = days,
                              freq = freq, numFreq = numFreq, delta = delta, 
                              startYear = startYear, endYear = endYear,
                              intraAnnualPeriod = c("wetSeason", "drySeason"), 
                              interAnnualPeriod = interAnnualPeriod, 
                              adhocPeriod = adhocPeriod)
  } else {
    # this is for matrix object
    message("Calculating output for class matrix...")
    days <- as.numeric(1:(ncol(data)-2))
    output <- get.sta.matrix(data = data, days = days,
                             freq = freq, numFreq = numFreq, delta = delta,
                             # startYear = startYear, endYear = endYear,
                             intraAnnualPeriod = c("wetSeason", "drySeason"),
                             interAnnualPeriod = interAnnualPeriod,
                             adhocPeriod = adhocPeriod, 
                             numCores = numCores, dirToSaveSTA = dirToSaveSTA)
    
    if( save ){
      if( is.null(dirToSaveSTA) ){
        message("'dirToSaveSTA' was not provided")
      } else {
        message(paste0("Saving output for class matrix. See directory: ",
                       dirToSaveSTA))
        
        mean <- output$sta$mean
        mean_basicStats <- output$sta$mean_basicStats
        annual <- output$sta$annual
        annual_basicStats <- output$sta$annual_basicStats
        semiannual <- output$sta$semiannual
        semiannual_basicStats <- output$sta$semiannual_basicStats
        
        save(mean, mean_basicStats, annual, annual_basicStats, 
             semiannual, semiannual_basicStats, 
             file = paste0( dirToSaveSTA, "/sta_matrix_output.RData" ))
      }
    } 
  }
  
  # sta.out <- structure(list(output = output), class = "sta")
  if(class(data) == "numeric"){
    sta.out <- structure(output, class = "staNumeric")
  } else {
    sta.out <- structure(output, class = "staMatrix")
  }
  
  if( plot ){
    message("Calculating plots...")
    
    if( class(sta.out) == "staNumeric"){
      plot.staNumeric(x = sta.out, startYear = startYear, endYear = endYear,
                      interAnnualPeriod = interAnnualPeriod, significance = significance)
    }
    
    if( class(sta.out) == "staMatrix" ){
      plot.staMatrix(x = sta.out, significance = significance, master = master)
    }
    
  } else {
    sta.out
  }
  
  # if( plot ){
  #   message("Calculating plots...")
  #   plot_sta(output = sta.out, startYear = startYear, endYear = endYear,
  #            interAnnualPeriod = interAnnualPeriod, significance = significance, master = master)
  #   # getPlot.sta.matrix(output = output, significance = significance, master = master)
  #   if( class(sta.out) == "staMatrix" & save == T ){
  #     # if( missing(master) ){
  #     #   master <- getMaster()
  #     # }
  #     message("Please wait, the 30 raster files associated with this STA are being saved")
  #     saveMaps(output = sta.out, master = master, dirToSaveSTA = dirToSaveSTA)
  #   }
  #   invisible(sta.out)
  # } else {
  #   sta.out
  # }
}
#
#' Plot method for \code{sta} function
#' 
#' This function returns a plot
#' 
#' @param x                 an object of class "staNumeric"
#' @param ...               additional plot parameters
#' @param startYear         numeric, time series initial year                            
#' @param endYear           numeric, time series last year                            
#' @param interAnnualPeriod numeric vector indicating the number of years to be considered in STA. For instance,
#'                          in a 19 years time series, 1:5, indicates that the first five years will be
#'                          utilized for STA. Similarly, c(2,6,10) indicates that the second, sith and tenth
#'                          years will be utilized for STA. See \code{\link[sta]{sta}}
#' @param significance      numeric indicating significance of each shape parameter trend      
#' 
#' @rdname plot.staNumeric
#' @method plot staNumeric
#' @export
#' 
#' @seealso \code{\link[sta]{sta}}
plot.staNumeric <- function(x, ...,  startYear, endYear, interAnnualPeriod,
                                significance){
  
  years <- startYear:endYear
  
  nf <- layout( matrix( c(1, 2, 3, 4), 4, 1, byrow = T ),
                c(12.5, 12.5, 12.5, 12.5), c(3, 3, 3, 3), T )
  
  # TEMP <- sta_wetSeason_21016$mean$globalDaysToAnalyze
  begin <- seq(1, length(x$intervalsUsedBasicStats), 2)
  end <- seq(2, length(x$intervalsUsedBasicStats), 2)
  intervals <- c(sapply(1:length(begin), function(s) 
    x$intervalsUsedBasicStats[begin[s]]:x$intervalsUsedBasicStats[end[s]]))
  
  daysForFit <- 1:length(x$data)
  fit <- rep(NA, length(daysForFit))
  fit[x$daysUsedFit] <- x$fit
  basicAnalysis <- fit
  basicAnalysis[-c(intervals)] <- NA
  
  ts_output <- ts(x$data, start = c(startYear, 1), end = c(endYear-1, x$freq), 
                  frequency = x$freq)
  fit_output <- ts(fit, start = c(startYear, 1), end = c(endYear-1, x$freq), 
                   frequency = x$freq)
  basic_output <- ts(basicAnalysis, start = c(startYear, 1), end = c(endYear-1, x$freq), 
                     frequency = x$freq)
  
  yRan <- range(c(x$data, fit, basicAnalysis), na.rm = T)
  
  draw_timeSeries <- c(4.1, 4.5, 1, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  par(mar = draw_timeSeries, cex.axis = 1.5, cex.lab= 1.5, font.lab = 2)
  plot(ts_output, ylim = yRan,
       type = "p", ylab = "NDMI", xlab = "", 
       col = "lightgray", lwd = 4)
  # cex.axis = 5, cex.lab = 5)
  par(new=T)
  plot(fit_output, ylim = yRan, lwd = 2,
       type = "l", ylab = "NDMI", xlab = "")
  # cex.axis = 5, cex.lab = 5)
  par(new=T)
  plot(basic_output, ylim = yRan,
       type = "p", ylab = "NDMI", xlab = "",
       col = "red", pch = 16)
  # font.lab = 2,
  # cex.axis = 5, cex.lab = 5)
  
  # --- MEAN
  
  COLOR <- "black"
  if(!is.null(significance)){
    if(x$sta$mean$pval < significance){
      COLOR <- "red"
    }
  }
  
  yRan <- range(c(x$sta$mean$harmCoeffs, 
                  x$sta$mean$linearTrend), na.rm = T)
  yRan[1] <- yRan[1] - 0.02
  yRan[2] <- yRan[2] + 0.02
  # draw_params_trends <- c(5.1, 4.1, 0, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  draw_params_trends <- c(4.1, 4.1, 1, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  par(mar = draw_params_trends)
  plot(x = 1:length(x$sta$mean$harmCoeffs), 
       x$sta$mean$harmCoeffs, type = "p", xlab = "", ylab = "mean",
       xaxt = "n",
       axes = F,
       ylim = yRan, pch = 16, cex = 2, cex.axis = 1.5, cex.lab = 1.5)
  axis(1, labels = rep("", length(x$sta$mean$harmCoeffs)),
       at = 1:length(x$sta$mean$harmCoeffs), cex.axis = 1.5)
  RAN <- range((x$sta$mean$harmCoeffs))
  axis(2, labels = round(seq(RAN[1], RAN[2], length = 4),2),
       at = seq(RAN[1], RAN[2], length = 4), cex.axis = 1.5)
  lines(x$sta$mean$linearTrend, lwd = 2)
  title( paste0("slope: ", round(x$sta$mean$slope, 3),
                "     ", "p-value: ", round(x$sta$mean$pval, 3)),
         col.main = COLOR )
  
  # --- ANNUAL
  
  COLOR <- "black"
  if(!is.null(significance)){
    if(x$sta$annual$pval < significance){
      COLOR <- "red"
    }
  }
  
  yRan <- range(c(x$sta$annual$harmCoeffs, 
                  x$sta$annual$linearTrend), na.rm = T)
  yRan[1] <- yRan[1] - 0.02
  yRan[2] <- yRan[2] + 0.02
  # draw_params_trends <- c(5.1, 4.1, 0, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  draw_params_trends <- c(4.1, 4.1, 1, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  par(mar = draw_params_trends)
  plot(x = 1:length(x$sta$annual$harmCoeffs), 
       x$sta$annual$harmCoeffs, type = "p", 
       xlab = "", ylab = "annual",
       xaxt = "n", axes = F,
       ylim = yRan, pch = 16, cex = 2, cex.axis = 1.5, cex.lab = 1.5)
  axis(1, labels = rep("", length(x$sta$annual$harmCoeffs)),
       at = 1:length(x$sta$annual$harmCoeffs), cex.axis = 1.5)
  RAN <- range(x$sta$annual$harmCoeffs)
  axis(2, labels = round(seq(RAN[1], RAN[2], length = 4),2),
       at = seq(RAN[1], RAN[2], length = 4), cex.axis = 1.5)
  lines(x$sta$annual$linearTrend, lwd = 2)
  title( paste0("slope: ", round(x$sta$annual$slope, 3), 
                "     ", "p-value: ", round(x$sta$annual$pval, 3)),
         col.main = COLOR)
  
  # --- SEMI-ANNUAL
  
  COLOR <- "black"
  if(!is.null(significance)){
    if(x$sta$semiannual$pval < significance){
      COLOR <- "red"
    }
  }
  
  yRan <- range(c(x$sta$semiannual$harmCoeffs, 
                  x$sta$semiannual$linearTrend), na.rm = T)
  yRan[1] <- yRan[1] - 0.02
  yRan[2] <- yRan[2] + 0.02
  # draw_params_trends <- c(5.1, 4.1, 0, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  draw_params_trends <- c(4.1, 4.1, 1, 2.1)  #  c(5.1, 5.1, 2.1, 2.1)
  par(mar = draw_params_trends)
  plot(x = 1:length(x$sta$semiannual$harmCoeffs), 
       x$sta$semiannual$harmCoeffs,
       type = "p", xlab = "", ylab = "semi-annual",
       xaxt = "n",
       axes = F,
       ylim = yRan, pch = 16, cex = 2, cex.axis = 1.5, cex.lab = 1.5)
  axis(1, labels = years[interAnnualPeriod],
       at = 1:length(x$sta$semiannual$harmCoeffs), cex.axis = 1.5)
  RAN <- range(x$sta$semiannual$harmCoeffs)
  axis(2, labels = round(seq(RAN[1], RAN[2], length = 4),2),
       at = seq(RAN[1], RAN[2], length = 4), cex.axis = 1.5)
  lines(x$sta$semiannual$linearTrend, lwd = 2)
  title( paste0("slope: ", round(x$sta$semiannual$slope, 3),
                "     ", "p-value: ", round(x$sta$semiannual$pval, 3)),
         col.main = COLOR )
}
#
#' Plot method for \code{sta} function
#' 
#' This function returns a \code{\link[mapview]{mapview}} object
#' 
#' @param x                 an object of class "staMatrix"
#' @param ...               additional plot parameters
#' @param significance      numeric indicating significance of each shape parameter trend
#' @param master            character with full path to a \emph{.tif} file which is used to transfer STA
#'                          output to a raster file.
#' 
#' @rdname plot.staMatrix
#' @method plot staMatrix
#' @export
#' 
#' @seealso \code{\link[sta]{sta}}
plot.staMatrix <- function(x, ..., significance, master){
  
  m_mean <- matrixToMapview(mat = x$sta$mean, significance = significance,
                            master = master,
                            color = "Reds", label = "mean")
  m_annual <- matrixToMapview(mat = x$sta$annual, significance = significance,
                              master = master,
                              color = "Greens", label = "annual")
  m_semiannual <- matrixToMapview(mat = x$sta$semiannual, 
                                  significance = significance,
                                  master = master,
                                  color = "Blues", label = "semiAnnual")
  
  maps <- m_mean + m_annual + m_semiannual
  
  maps@map
}




