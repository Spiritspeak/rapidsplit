#' @title rapidsplithalf package
#' @name rapidsplithalf
#' @aliases rapidsplithalf-package
#' @import Rcpp doParallel foreach
#' @useDynLib rapidsplithalf
#' @importFrom stats cor setNames qnorm pt pnorm quantile weighted.mean
#' @importFrom grDevices rgb
#' @importFrom graphics text
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom fastmatch fmatch %fin%
#' @importFrom kit funique countOccur uniqLen
#' @author Sercan Kahveci
#' @description
#' To learn more about rapidsplithalf, view the introductory vignette:
#' \code{vignette("rapidsplithalf",package="rapidsplithalf")}
#' 
NULL