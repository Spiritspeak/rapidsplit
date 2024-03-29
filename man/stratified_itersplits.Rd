% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{stratified_itersplits}
\alias{stratified_itersplits}
\title{stratified_itersplits}
\usage{
stratified_itersplits(itercount, groupsizes)
}
\arguments{
\item{itercount}{number of iterations}

\item{groupsizes}{vector of number of RTs per group to stratify}
}
\value{
A matrix with zeroes and ones
}
\description{
generate stratified splits for a single participant
}
\details{
This first equally splits what can be equally split within groups.
Then it randomly splits all the leftovers.
}
