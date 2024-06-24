#' @name rapidsplit
#' @aliases rapidsplit-package
#' @import Rcpp doParallel foreach
#' @useDynLib rapidsplit
#' @importFrom stats cor setNames qnorm pt pnorm quantile weighted.mean
#' @importFrom grDevices rgb
#' @importFrom graphics text
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
NULL