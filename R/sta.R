#' Statistical Seasonal Trend Analysis for numeric vector or \code{RasterStack}
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
#' @param significance      numeric in the interval (0,1) to assess statistical significance of trend analysis.
#'                          \code{NULL} by default.                           
#' @param save              logical, should STA output be saved, default is \code{FALSE}
#' @param dirToSaveSTA      character with full path name used to save \code{sta} progress report. When 
#'                          \code{save = TRUE}, \code{sta}'s output is saved on this path.
#' @param numCores          integer given the number of cores to use; pertinent when \code{data} is a
#'                          RasterStack or a \code{\link[base]{matrix}}
#'                          
#' @note STA is based on the following ideas. Let \eqn{y_t} denote the value of a periodic time
#' series at time-point \eqn{t}. It is well-known that this type of observations can be modeled
#' as:
#' 
#' \eqn{y_t = a_0 + a_1 cos( (2 \pi t)/L - \phi_1) + ... +  a_K cos( (2 \pi K t)/L - \phi_K) + \varepsilon_t}, \eqn{t=1,\ldots,L}.
#'                 
#' This model is known as harmonic regression. \eqn{L} denotes the number of observations per period, \eqn{K} is the number of 
#' harmonics included in the fit, \eqn{a_i}'s and \eqn{\phi_i}'s are called amplitude coefficients and phase angles,
#' respectively. \eqn{K} can be known empirically. Amplitudes and phase angle can be obtained as the solution of a least-squares
#' problem. 
#' 
#' In vegetation monitoring, amplitudes and phase angles are known as \emph{shape parameters}. In particular,
#' \eqn{a_0}, \eqn{a_1} and \eqn{a_2} are called \emph{mean} and \emph{annual} and \emph{semiannual} amplitudes, respectively.
#' Applying the harmonic regression model to observations over \eqn{P} periods of length \eqn{L} each, results
#' in estimates of shape parameters for each period. Thus, focusing only on amplitudes, \code{\link{sta}} yields
#' time series of mean, annual and semiannual parameters. A subsequent Mann-Kendall test for trend is performed
#' on each of these series.
#' 
#' @export
#' 
#' @examples
#' sta_marismas <- sta(data=marismas, freq=36)
#' str(sta_marismas)
#' plot(sta_marismas)
#' plot(sta_marismas, significance=0.09)
#' 
#' # Use of interAnnualPeriod
#' sta_21016 <- sta(data = marismas, freq = 36, interAnnualPeriod = c(2, 10, 16))
#' plot(sta_21016)
#' 
#' # Use of intraAnnualPeriod
#' sta_drySeason_218 <- sta(data = marismas, freq = 36,
#'                      interAnnualPeriod = 2:18, intraAnnualPeriod = "drySeason")
#' plot(sta_drySeason_218)
#' 
#' # Use of adhocPeriod and significance
#' adhoc <- list()
#' beginPeriod <- (1:17) * 36
#' endPeriod <- 2:18 * 36 
#' adhoc$partial <- c( sapply(1:length(beginPeriod), 
#'                  function(s) c(beginPeriod[s]+1, endPeriod[s]) ) )
#' adhoc$full <- c( sapply(1:length(beginPeriod), 
#'               function(s) (beginPeriod[s]+1):endPeriod[s]) )
#' sta_adhoc_218 <- sta(data = marismas, freq = 36, interAnnualPeriod = 2:18,
#'                  startYear = 2000, endYear = 2018, adhocPeriod = adhoc, significance=0.05)
#' plot(sta_adhoc_218)
#' 
#' # Use of ndmi RasterStack
#' \donttest{
#' ndmi_path = system.file("extdata", "ndmi.tif", package = "sta")
#' ndmiSTACK <- stack(ndmi_path)
#' dir.create(path=paste0(system.file("extdata", package="sta"), "/output_ndmi"),
#'           showWarnings=FALSE)
#' outputDIR = paste0(system.file("extdata", package="sta"), "/output_ndmi")
#' 
#' sta_ndmi_21016 <- sta(data = ndmiSTACK, freq = 36,
#'                   numFreq = 4, delta = 0.2, intraAnnualPeriod = "wetSeason",
#'                   startYear = 2000, endYear = 2018, interAnnualPeriod = c(2,10,16),
#'                   save = TRUE, numCores = 2L, dirToSaveSTA = outputDIR)}
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
#' to a \code{RasterStack} has this format. 
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
#' or \code{\link[base]{matrix}}, \code{interAnnualPeriod = 1:((ncol(data)-2)/freq)}.
#' 
#' Since \code{adhocPeriod} defines an inter annual period "ad-hoc", the specific days of this ad-hoc
#' season must be known in advance and consequently, the specific time-points (with respect to the 
#' time series under consideration) must be provided in a numeric vector.
#' 
#' When \code{save=T}, a valid \code{dirToSaveSTA} must be provided, that is, this folder should have been
#' created previously. In this case, \code{sta}'s output is saved on \code{dirToSaveSTA}. This version
#' saves arrays of STA of the mean, annual and semi-annual parameters (along with their corresponding basic statistics)
#' in the file \code{sta_matrix_output.RData} inside \code{dirToSaveSTA}. Also, in the same directory,
#' the file \code{sta_progress.txt} records the progress of the STA process.
#' 
#' \code{save=T}, \code{dirToSaveSTA}, \code{numCores} and \code{master} are required when \code{data} is either a
#' \code{RasterStack} or a \code{\link[base]{matrix}}. The aforementioned basic statistics are: mean and standard deviation
#' of the time series of annual maximum and minimum as well as the global minima and maxima.
#' 
#' @return When \code{class(data)} is a \code{numeric}, an object of class "staNumeric" containing:
#' \item{data}{numeric vector}
#' \item{freq}{integer with the number of observations per period}
#' \item{startYear}{numeric, time series initial year}
#' \item{endYear}{numeric, time series last year}
#' \item{intraAnnualPeriod}{character indicating seasons (wet or dry)}
#' \item{interAnnualPeriod}{numeric vector indicating the number of years considered in STA}
#' \item{ts}{time series object; \code{data} in \code{\link[stats]{ts}} format}
#' \item{fit}{numeric vector with output of \code{\link[geoTS]{haRmonics}}}
#' \item{sta}{a list containing the following elements:}
#' \itemize{
#' \item \code{mean} a list of 12 elements with STA output for shape parameter \emph{mean}
#' \item \code{annual} a list of 12 elements with STA output for shape parameter \emph{annual}
#' \item \code{semiannual} a list of 12 elements with STA output for shape parameter \emph{semiannual}
#' }
#' \item{significance}{numeric in the interval (0,1) or \code{NULL} when default used}
#'
#' When \code{class(data)} is a \code{RasterStack} or a \code{\link[base]{matrix}}, an object of class
#' "staMatrix" containing:
#' \item{freq}{integer with the number of observations per period}
#' \item{daysUsedFit}{integer vector indicating the indices used to fit harmonic regression. see \code{\link[geoTS]{haRmonics}}}
#' \item{intervalsUsedBasicStats}{numeric vector indicating the indices used on calculation of basic statistics}
#' \item{sta}{a list containg the following elements:}
#' \itemize{
#' \item \code{mean} a matrix of 4 columns with STA output for shape parameter \emph{mean}. First two columns provide geolocation of analyzed pixels, 
#' third and fourth columns show p-value and slope estimate of STA, respectively
#' \item \code{mean_basicStats} a matrix of 10 columns with basic statistics for shape parameter \emph{mean}. First two columns provide geolocation of analyzed pixels, 
#' from third to tenth columns show mean, standard deviation, global minimum, and maximum of minimum values, as well as mean, standard deviation,
#' global minimum, and maximum of maximum values, respectively
#' \item \code{annual} a matrix of 4 columns with STA output for shape parameter \emph{annual}. First two columns provide geolocation of analyzed pixels, 
#' third and fourth columns show p-value and slope estimate of STA, respectively
#' \item \code{annual_basicStats} a matrix of 10 columns with basic statistics for shape parameter \emph{annual}. First two columns provide geolocation of analyzed pixels, 
#' from third to tenth columns show mean, standard deviation, global minimum, and maximum of minimum values, as well as mean, standard deviation,
#' global minimum, and maximum of maximum values, respectively
#' \item \code{semiannual} a matrix of 4 columns with STA output for shape parameter \emph{semiannual}. First two columns provide geolocation of analyzed pixels, 
#' third and fourth columns show p-value and slope estimate of STA, respectively
#' \item \code{semiannual_basicStats} a matrix of 10 columns with basic statistics for shape parameter \emph{semiannual}. First two columns provide geolocation of analyzed pixels, 
#' from third to tenth columns show mean, standard deviation, global minimum, and maximum of minimum values, as well as mean, standard deviation,
#' global minimum, and maximum of maximum values, respectively 
#' }
#' @references Eastman, R., Sangermano, F., Ghimine, B., Zhu, H., Chen, H., Neeti, N., Cai, Y., Machado E., Crema, S. (2009).
#' \emph{Seasonal trend analysis of image time series},
#' International Journal of Remote Sensing, \bold{30(10)}, 2721--2726.
#' 
sta <- function(data, freq, numFreq = 4, delta = 0, startYear = 2000, endYear = 2018,
                intraAnnualPeriod = c("wetSeason", "drySeason"), interAnnualPeriod, adhocPeriod = NULL,
                significance = NULL, save = FALSE, dirToSaveSTA = NULL, numCores = 20){

  if( !( inherits(data, "numeric") | inherits(data, "RasterStack") | inherits(data, "matrix") ) ){
    stop("'data' must be either numeric, RasterStack or matrix")
  }
  
  if( inherits(data, "RasterStack") ){ # class(data) == "RasterStack"
    message("Transferring 'data' to matrix...")
    data <- rasterToPoints(data)
  }
  
  if( missing(freq) ){
    stop("'freq' must be provided")
  }
  
  intraAnnualPeriod <- match.arg(intraAnnualPeriod)
  
  if( missing(interAnnualPeriod) ){
    message("Calculating 'interAnnualPeriod'...")
    if( inherits(data, "numeric") ){ #  class(data) == "numeric"
      interAnnualPeriod <- 1:(length(data)/freq)
    } else {
      interAnnualPeriod <- 1:((ncol(data)-2)/freq)
    }
  } else {
    if( inherits(data, "numeric") ){ #  class(data) == "numeric"
      if( length(interAnnualPeriod) > length(data)/freq )
       stop("length(data)/freq must be greater or equal than length(interAnnualPeriod)")
    } else {
      if( length(interAnnualPeriod) > (ncol(data)-2)/freq )
        stop("(ncol(data)-2)/freq must be greater or equal than length(interAnnualPeriod)")
    }
  }
  
  if( !is.null(adhocPeriod) ){
    if( !is.list(adhocPeriod) | length(adhocPeriod) != 2 ){
      stop("'adhocPeriod' must be a list of lenght 2")
    }
  }
  
  if( inherits(data, "numeric") ){ #class(data) == "numeric"
    message("Calculating output for class numeric...")
    output <- get.sta.numeric(data = data, #days = days,
                              days = as.numeric(1:length(data)),
                              freq = freq, numFreq = numFreq, delta = delta, 
                              startYear = startYear, endYear = endYear,
                              intraAnnualPeriod = intraAnnualPeriod, 
                              interAnnualPeriod = interAnnualPeriod, 
                              adhocPeriod = adhocPeriod)
    output$significance <- significance
  } else {
    # this is for matrix object
    message("Calculating output for class matrix...")
    output <- get.sta.matrix(data = data, 
                             days = as.numeric(1:(ncol(data)-2)),
                             freq = freq, numFreq = numFreq, delta = delta,
                             startYear = startYear, endYear = endYear,
                             intraAnnualPeriod = intraAnnualPeriod, 
                             interAnnualPeriod = interAnnualPeriod,
                             adhocPeriod = adhocPeriod, 
                             numCores = numCores, dirToSaveSTA = dirToSaveSTA)
    output$significance <- significance
    
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
  
  if( inherits(data, "numeric") ){# class(data) == "numeric"
    sta.out <- structure(output, class = "staNumeric")
  } else {
    sta.out <- structure(output, class = "staMatrix")
  }
  
  sta.out
}

#' Plot method for \code{sta} function
#' 
#' This function returns a plot
#' 
#' @param x                 an object of class "staNumeric"
#' @param significance      numeric indicating significance of each shape parameter trend   
#' @param ...               additional plot parameters    
#' 
#' @rdname plot.staNumeric
#' @method plot staNumeric
#' @export
#' 
#' @seealso \code{\link{sta}}
plot.staNumeric <- function(x, 
                            # startYear, endYear, interAnnualPeriod=1:(length(startYear:endYear)),
                            significance=NULL, ...){
  
  if( !inherits(x, "staNumeric") ){ # class(x) != "staNumeric"
    stop("'x' must be of class staNumeric")
  }
  
  if(!is.null(significance)){
    x$significance <- significance
  }
  
  years <- x$startYear:x$endYear
  
  nf <- layout( matrix( c(1, 2, 3, 4), 4, 1, byrow = T ),
                c(12.5, 12.5, 12.5, 12.5), c(3, 3, 3, 3), T )
  
  begin <- seq(1, length(x$intervalsUsedBasicStats), 2)
  end <- seq(2, length(x$intervalsUsedBasicStats), 2)
  intervals <- c(sapply(1:length(begin), function(s) 
    x$intervalsUsedBasicStats[begin[s]]:x$intervalsUsedBasicStats[end[s]]))
  
  daysForFit <- 1:length(x$data)
  fit <- rep(NA, length(daysForFit))
  fit[x$daysUsedFit] <- x$fit
  basicAnalysis <- fit
  basicAnalysis[-c(intervals)] <- NA
  
  ts_output <- ts(x$data, start = c(x$startYear, 1), end = c(x$endYear, x$freq), 
                  frequency = x$freq)
  fit_output <- ts(fit, start = c(x$startYear, 1), end = c(x$endYear, x$freq), 
                   frequency = x$freq)
  basic_output <- ts(basicAnalysis, start = c(x$startYear, 1), end = c(x$endYear, x$freq), 
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

  # --- MEAN
  
  COLOR <- "black"
  if(!is.null(x$significance)){
    if(x$sta$mean$pval < x$significance){
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
  if(!is.null(x$significance)){
    if(x$sta$annual$pval < x$significance){
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
  if(!is.null(x$significance)){
    if(x$sta$semiannual$pval < x$significance){
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
  # axis(1, labels = years[interAnnualPeriod],
  #      at = 1:length(x$sta$semiannual$harmCoeffs), cex.axis = 1.5)
  axis(1, labels = years[x$interAnnualPeriod],
       at = 1:length(x$interAnnualPeriod), cex.axis = 1.5)
  # axis(1, labels = rep("", length(x$sta$semiannual$harmCoeffs)),
  #      at = 1:length(x$sta$semiannual$harmCoeffs), cex.axis = 1.5)
  RAN <- range(x$sta$semiannual$harmCoeffs)
  axis(2, labels = round(seq(RAN[1], RAN[2], length = 4),2),
       at = seq(RAN[1], RAN[2], length = 4), cex.axis = 1.5)
  lines(x$sta$semiannual$linearTrend, lwd = 2)
  title( paste0("slope: ", round(x$sta$semiannual$slope, 3),
                "     ", "p-value: ", round(x$sta$semiannual$pval, 3)),
         col.main = COLOR )
}

#' Plot method for \code{sta} function
#' 
#' This function displays some maps of \code{\link[mapview]{mapview-class}}
#' 
#' @param x                 an object of class "staMatrix"
#' @param significance      numeric indicating significance of each shape parameter trend
#' @param master            RasterLayer used to transfer STA output to raster layers.
#' @param ...               additional plot parameters
#' 
#' @rdname plot.staMatrix
#' @method plot staMatrix
#' @export
#' 
#' @seealso \code{\link[sta]{sta}}, \code{\link[sta]{getMaster}}, \code{\link[geoTS]{matrixToRaster}}
plot.staMatrix <- function(x, significance=NULL, master, ...){
  if( !inherits(x, "staMatrix") ){ #class(x) != "staMatrix"
    stop("'x' must be of class staMatrix")
  }
  
  if(!is.null(significance)){
    x$significance <- significance
  }
  
  if(missing(master)){
    stop("'master' must be provided")
  }
  
  m_mean <- matrixToMapview(mat = x$sta$mean, significance = x$significance,
                            master = master,
                            color = "Reds", label = "mean")
  m_annual <- matrixToMapview(mat = x$sta$annual, significance = x$significance,
                              master = master,
                              color = "Greens", label = "annual")
  m_semiannual <- matrixToMapview(mat = x$sta$semiannual, 
                                  significance = x$significance,
                                  master = master,
                                  color = "Blues", label = "semiAnnual")
  
  maps <- m_mean + m_annual + m_semiannual
  
  maps@map
}




