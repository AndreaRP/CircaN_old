#' circan
#'
#' Nonlinear least squares model for accurate detection of circadian expression patterns.
#' @param data Dataframe containing the expression data. Samples must be in columns and genes in rows.
#' For an example see data(expression_example).
#' @param s2c Dataframe containing the metadata for the samples. Must have at least a 'sample' column
#' with the sample name as it appears in the data matrix;
#' a 'time' column with the time point the sample was collected;
#' and an 'ind' column containing information for the individual the sample comes from.
#' For an example see data(metadata_example).
#' @param shiny Is the package running in a shiny app? default to FALSE.
#' @param mode Algorithm to use in the NLS regression. Must be one of 'default' for Gauss-Newton, 'plinear' for the Golub-Pereyra algorithm
#' for partially linear least-squares models and 'port' for the ‘nl2sol’ algorithm from the Port library. Default is default. See nls documentation
#' for extended info.
#' @param init_value Initial value for the period. Default is set to 24.
#' @param max_per Maximum period to regress. Default is set to Inf.
#' @param min_per Minimum period to regress. Default is set to -Inf.
#' @keywords CircaN circadian regression
#' @export
#' @examples
#' # This runs CircaN on the example data with the 'Port' algorithm.
#' circan(data=expression_example, s2c=metadata_example, mode="port")

circan <- function(data, s2c, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
  s2c$time <- as.numeric(as.character(s2c$time))
  s2c$ind <- as.integer(as.character(s2c$ind))
  s2c$sample <- as.character(s2c$sample)
  # metadata must have sample, ind and time columns. Data is a data frame or matrix
  # with the ids in the first column
  t <- unique(as.numeric(as.character(s2c$time)))
  data0 <- data
  total_genes <- nrow(data)
  rownames(data) <- as.character(data[,1])
  data <- data[,-1]
  data <- as.matrix(data)
  results.cols <- c("feature"
                    ,"estimate.amp","std.error.amp","statistic.amp","p.value.amp"
                    ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
                    ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
                    , "r")
  results <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), results.cols)
  for(gene in 1:nrow(data)){
    # Standarize data
    data.st <- (data[gene,]-mean(data[gene,]))/sqrt(var(data[gene,]))
    # Group data by timepoint
    df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
    df0 <- dplyr::rename(df0, data = x)
    df <- as.data.frame(df0[,c("data", "ind", "time")])
    df$data <- as.numeric(df$data)
    df$time <- as.numeric(df$time)
    gd <- nlme::groupedData(data~time|ind,data=df)
    # plot(gd)
    result = tryCatch({
      # Try fitting all parameters
      nls.model=nls(data~amp*sin(2*pi*time/per+phase)
                    , start=list(amp=1, phase=0, per=init_value)
                    , algorithm = mode
                    , lower = list(per=min_per)
                    , upper = list(per=max_per)
                    , data = gd
                    , control = list(maxiter = 200)
      )
      stats <- as.data.frame(broom::tidy(nls.model))
      rownames(stats) <- stats[,1]
      stats <- stats[,-1]
      vec <- as.numeric(c(t(stats)))
      names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
      # R (goodness of fit)
      r <- cor(gd$data, predict(nls.model))
    }, error = function(e1) {
      convergence_error1 <- grepl("Convergence failure", as.character(e1))
      if(convergence_error1){
        tryCatch({
          # Try with fixed period
          nls.model=nls(  data~amp*sin(2*pi*time/init_value+phase)
                          , start=list(amp=1,phase=0)
                          , algorithm = mode
                          , data=gd
                          , control = list(maxiter = 200
                                           # , warnOnly = TRUE
                          )
          )
          stats <- as.data.frame(broom::tidy(nls.model))
          rownames(stats) <- stats[,1]
          stats <- stats[,-1]
          vec <- as.numeric(c(t(stats)))
          names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
          # As we haven't estimated the period we add the data to vec
          temp <- c(24, "NA", "NA", "NA")
          names(temp) <- c("estimate.period", "std.error.period"
                           , "statistic.period", "p.value.period")
          vec <<- c(vec, temp)
          # R (goodness of fit)
          r <<- cor(gd$data, predict(nls.model))
        }, error = function(e2) {
          convergence_error2 <- grepl("Convergence failure", as.character(e2))
          if(convergence_error2){
            # Can't be adjusted to a circadian curve
            vec <<- rep("NA", times = ncol(results)-1)
            r <<- "NA"
          }
        })
      }
    }, finally = {
      res <- rbind(c(feature=rownames(data)[gene], vec, r = r))
      results[gene,] <- res
    })
    if(shiny) incProgress(amount=1/total_genes)
  }
  # Calc adjusted p.val for period
  results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method="BH")
  results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method="BH")
  results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method="BH")
  return(results)
}



# circan_old <- function(data, s2c, shiny=FALSE, mode="default", init_value=24, max_per=Inf, min_per=-Inf){
#   s2c$time <- as.numeric(as.character(s2c$time))
#   s2c$ind <- as.integer(as.character(s2c$ind))
#   s2c$sample <- as.character(s2c$sample)
#   # metadata must have sample, ind and time columns. Data is a data frame or matrix
#   # with the ids in the first column
#   t <- unique(as.numeric(as.character(s2c$time)))
#   data0 <- data
#   total_genes <- nrow(data)
#   rownames(data) <- as.character(data[,1])
#   data <- data[,-1]
#   data <- as.matrix(data)
#   results.cols <- c("feature"
#                     ,"estimate.amp","std.error.amp","statistic.amp","p.value.amp"
#                     ,"estimate.phase","std.error.phase","statistic.phase","p.value.phase"
#                     ,"estimate.per" ,"std.error.per","statistic.per", "p.value.per"
#                     , "r")
#   results <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), results.cols)
#   for(gene in 1:nrow(data)){
#     # Standarize data
#     data.st <- (data[gene,]-mean(data[gene,]))/sqrt(var(data[gene,]))
#     # Group data by timepoint
#     df0 <- merge(data.st, s2c, by.x = 0, by.y = "sample")
#     df0 <- dplyr::rename(df0, data = x)
#     df <- as.data.frame(df0[,c("data", "ind", "time")])
#     df$data <- as.numeric(df$data)
#     df$time <- as.numeric(df$time)
#     gd <- nlme::groupedData(data~time|ind,data=df)
#     # plot(gd)
#     result = tryCatch({
#       # Fit
#       nls.model=nls(data~amp*sin(2*pi*time/per+phase)
#                     , start=list(amp=1, phase=0, per=init_value)
#                     , algorithm = mode
#                     , lower = list(per=min_per)
#                     , upper = list(per=max_per)
#                     , data = gd
#                     , control = list(maxiter = 200)
#       )
#       stats <- as.data.frame(broom::tidy(nls.model))
#       rownames(stats) <- stats[,1]
#       stats <- stats[,-1]
#       vec <- as.numeric(c(t(stats)))
#       names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#       # R (goodness of fit)
#       r <- cor(gd$data, predict(nls.model))
#     }, error = function(e) {
#       # print(e)
#       result2 = tryCatch({
#         # Try with fixed period
#         nls.model=nls(  data~amp*sin(2*pi*time/init_value+phase)
#                         , start=list(amp=1,phase=0)
#                         , algorithm = mode
#                         , data=gd
#                         , control = list(maxiter = 200
#                                          # , warnOnly = TRUE
#                         )
#         )
#         stats <- as.data.frame(broom::tidy(nls.model))
#         rownames(stats) <- stats[,1]
#         stats <- stats[,-1]
#         vec <- as.numeric(c(t(stats)))
#         names(vec) <- c(outer(colnames(stats), rownames(stats), paste, sep = "."))
#         # As we haven't estimated the period we add the data to vec
#         temp <- c(24, "NA", "NA", "NA")
#         names(temp) <- c("estimate.period", "std.error.period"
#                          , "statistic.period", "p.value.period")
#         vec <<- c(vec, temp)
#         # R (goodness of fit)
#         r <<- cor(gd$data, predict(nls.model))
#       }, error = function(e) {
#         # Can't be adjusted to a circadian curve
#         vec <<- rep("NA", times = ncol(results)-1)
#         r <<- "NA"
#       })
#     }, finally = {
#       res <- rbind(c(feature=rownames(data)[gene], vec, r = r))
#       results[gene,] <- res
#     })
#     if(shiny) incProgress(amount=1/total_genes)
#   }
#   # Calc adjusted p.val for period
#   results$BH.q.value.per <- p.adjust(as.numeric(results$p.value.per), method="BH")
#   results$BH.q.value.amp <- p.adjust(as.numeric(results$p.value.amp), method="BH")
#   results$BH.q.value.phase <- p.adjust(as.numeric(results$p.value.phase), method="BH")
#   return(results)
# }
