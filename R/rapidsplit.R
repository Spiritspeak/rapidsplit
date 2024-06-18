
cols2ids<-function(df){
  for(col in seq_along(df)){
    df[,col]<-as.numeric(as.factor(df[[col]]))
  }
  df
}

clamp.range<-function(x){
  x[x==min(x)]<- -1
  x[x==max(x)]<- 1
  x
}

makeCustomAggregator<-function(fun){
  function(x,mask){
    out<-numeric(ncol(mask))
    for(i in seq_len(ncol(mask))){
      out[i]<-fun(x[mask[,i,drop=T]])
    }
    return(out)
  }
}
# check for zero variances when standardizing
# accept only 1 subscore var, check for input lengths
datachecker<-function(data,subjvar,diffvars,stratvars,subscorevar,aggvar,standardize){
  if(!is.data.frame(data)){
    stop("Input dataset not a data.frame.")
  }
  if(length(subjvar)!=1 | length(aggvar)!=1 | !is.character(subjvar) | !is.character(aggvar)){
    stop("subjvar and aggvar must be character vectors of length 1.")
  }
  if(length(subscorevar)>1){
    stop("Can only specify 1 subscorevar.")
  }else if(length(subscorevar==1)){
    if(!is.character(subscorevar)){
      stop("subscorevar must be a character vector of length 1.")
    }
  }
  if(length(diffvars)>0){
    if(!is.character(diffvars)){
      stop("diffvars must be a character vector.")
    }
  }
  if(length(stratvars)>0){
    if(!is.character(stratvars)){
      stop("stratvars must be a character vector.")
    }
  }
  
  if(!all(c(subjvar,diffvars,stratvars,subscorevar,aggvar) %in% names(data))){
    stop("Not all specified columns present in data.")
  }
  if(length(intersect(diffvars,stratvars))>0){
    stop("Cannot have a variable specified as both diffvars and stratvars: ",
         intersect(diffvars,stratvars))
  }
  if(length(intersect(subscorevar,stratvars))>0){
    stop("Cannot have a variable specified as both subscorevar and stratvars: ",
         intersect(subscorevar,stratvars))
  }
  if(length(intersect(diffvars,subscorevar))>0){
    stop("Cannot have a variable specified as both diffvars and subscorevar: ",
         intersect(diffvars,subscorevar))
  }
  if(anyNA(data[[aggvar]])){
    stop(sum(is.na(data[[aggvar]]))," missing values present in variable to be aggregated.",
         " Please remove them first.")
  }
  if(!is.numeric(data[[aggvar]])){
    stop("Aggregation variable is not numeric.")
  }
  nvalues<-sapply(data[,diffvars,drop=FALSE],function(x){ length(unique(x))})
  if(any(nvalues!=2)){
    stop("Less or more than 2 unique values present in diffvars ",diffvars[nvalues!=2],
         ". The variables specified by this argument should feature two unique values, ",
         "indicating which condition (the lower value) ",
         "will be subtracted from the other (the higher value).")
  }
  
  if(length(subscorevar)>0){
    indices<-as.data.frame(table(unique(data[,c(subjvar,subscorevar,diffvars)])))
    subscoresubvalues<-table(do.call(paste,indices[,c(subjvar,subscorevar),drop=FALSE]))
    if(!all(subscoresubvalues == 2^length(diffvars))){
      stop("Some participants do not have data for all specified conditions within subscores: ",
           paste0("participant ",indices[[subjvar]][subscoresubvalues!=2^length(diffvars)],
                 " within subscore ",indices[,subscorevar][subscoresubvalues!=2^length(diffvars)],
                 collapse=", "))
    }
    condpaste<-do.call(paste,data[,c(subjvar,subscorevar,diffvars),drop=FALSE])
    condcounts<-table(condpaste)
    smallconds<-names(condcounts[condcounts<2])
    if(!all(condcounts>=2)){
      stop("Insufficient data (<2 obs) in 1 or more conditions belonging to these participants: ",
           paste0(unique(paste0(data[[subjvar]][condpaste %in% smallconds],
                               " within subscore ",
                               data[[subscorevar]][condpaste %in% smallconds])),
                 collapse=", "))
    }
    datapersubject<-table(do.call(paste,data[,c(subjvar,subscorevar),drop=FALSE]))
    if(standardize & !all(datapersubject>=4)){
      stop("Cannot compute some participants' standard deviation due to insufficient data: ",
           paste(names(datapersubject)[datapersubject<4],collapse=", "))
    }
    if(standardize){
      indices<-do.call(paste,data[,c(subjvar,subscorevar),drop=FALSE])
      nuniques<-tapply(X=data[[aggvar]],INDEX=indices,FUN=function(x){length(unique(x))})
      toofew<-names(nuniques[nuniques==1])
      if(any(nuniques==1)){
        stop("Cannot compute the SD of every participant/subscore, ",
             "due to less than 2 unique values in these participants/subscores: ",
             paste0(unique(paste0("participant ",data[[subjvar]][indices %in% toofew],
                                  " within subscore ",data[[subscorevar]][indices %in% toofew])),
                    collapse=", "))
      }
    }
  }else{
    indices<-as.data.frame(table(unique(data[,c(subjvar,diffvars)])))
    subscoresubvalues<-table(indices[[subjvar]])
    if(!all(subscoresubvalues == 2^length(diffvars))){
      stop("Some participants do not have data for all specified conditions: ",
           paste(unique(indices[[subjvar]][subscoresubvalues!=2^length(diffvars)]),collapse=", "))
    }
    condpaste<-do.call(paste,data[,c(subjvar,diffvars),drop=FALSE])
    condcounts<-table(condpaste)
    smallconds<-names(condcounts[condcounts<2])
    if(!all(condcounts>=2)){
      stop("Insufficient data (<2 obs) in 1 or more conditions belonging to these participants: ",
           paste(unique(data[[subjvar]][condpaste %in% smallconds]),collapse=", "))
    }
    datapersubject<-table(data[[subjvar]])
    if(standardize & !all(datapersubject>=4)){
      stop("Cannot compute some participants' standard deviation due to insufficient data: ",
           paste(names(datapersubject)[datapersubject<4],collapse=", "))
    }
    if(standardize){
      nuniques<-tapply(X=data[[aggvar]],INDEX=data[[subjvar]],FUN=function(x){length(unique(x))})
      toofew<-names(nuniques[nuniques==1])
      if(any(nuniques==1)){
        stop("Cannot compute the SD of every participant/subscore, ",
             "due to less than 2 unique values in these participants: ",
             paste0(toofew,collapse=", "))
      }
    }
  }
}

#' rapidsplit
#' @description A very fast algorithm for computing stratified permutated split-half reliability.
#'
#' @param data Dataset, a \code{data.frame}.
#' @param subjvar Subject ID variable name, a \code{character}.
#' @param diffvars Names of variables that determine which conditions 
#' need to be subtracted from each other, \code{character}.
#' @param stratvars Additional variables that the splits should 
#' be stratified by; a \code{character}.
#' @param subscorevar A \code{character} variable identifying subgroups within 
#' a participant's data from which separate scores should be computed. 
#' To compute a participant's final score, these subscores will be averaged together. 
#' A typical use case is the D-score of the implicit association task.
#' @param aggvar Name of variable whose values to aggregate, a \code{character}.
#' Examples include reaction times and error rates.
#' @param splits Number of split-halves to average, an \code{integer}.
#' @param aggfunc The function by which to aggregate the variable 
#' defined in \code{aggvar}; can be \code{"means"}, \code{"medians"}, 
#' or a \code{"custom"} function, given with \code{"custom.aggfunc"}.
#' @param custom.aggfunc A function that takes a numeric vector and outputs a single value,
#' to be used to aggregate a participant's values per condition. 
#' Only if \code{aggfunc} is set to \code{"custom"}.
#' @param standardize Whether to divide by scores by the subject's SD; a \code{logical}.
#' @param include.scores Include all individual split-half scores?
#' @param verbose Display progress bars?
#' @param check Check input for possible problems?
#'
#' @return A \code{list} containing 
#' 
#' * \code{r}: the averaged reliability.
#' 
#' * \code{allcors}: a vector with the reliability of each iteration.
#' 
#' * \code{nobs}: the number of participants.
#' 
#' * \code{scores}: the individual participants scores in each split-half, 
#' contained in a list with two matrices (Only included if requested with \code{include.scores}).
#' 
#' @export
#'
#' @examples 
#' 
#' data(foodAAT)
#' # Reliability of the double difference score:
#' # [RT(push food)-RT(pull food)] - [RT(push object)-RT(pull object)]
#' 
#' frel<-rapidsplit(data=foodAAT,
#'                  subjvar="subjectid",
#'                  diffvars=c("is_pull","is_target"),
#'                  stratvars="stimid",
#'                  aggvar="RT",
#'                  splits=1000)
#'                  
#' print(frel)
#' 
#' plot(frel,type="all")
#' 
#' # Unstratified reliability of the median RT
#' rapidsplit(data=foodAAT,
#'            subjvar="subjectid",
#'            aggvar="RT",
#'            splits=1000,
#'            aggfunc="medians")
#'            
#' # Compute a single random split-half reliability of the error rate
#' rapidsplit(data=foodAAT,
#'            subjvar="subjectid",
#'            aggvar="error",
#'            splits=1,
#'            aggfunc="means")
#'            
#' # Compute standardized mean differences between pulling and pushing RTs
#' # separately for both blocks. Then average the value from both blocks together.
#' rapidsplit(data=foodAAT,
#'            subjvar="subjectid",
#'            diffvars="is_pull",
#'            subscorevar="blocknum",
#'            aggvar="RT",
#'            splits=1000,
#'            standardize=TRUE)
#' 
rapidsplit<-function(data,subjvar,diffvars=NULL,stratvars=NULL,subscorevar=NULL,
                     aggvar,splits,
                     aggfunc=c("means","medians","custom"),custom.aggfunc=NULL,
                     standardize=FALSE,include.scores=TRUE,verbose=TRUE,check=TRUE){
  aggfunc<-match.arg(aggfunc)
  if(aggfunc=="means"){
    aggfunc<-meansByMask
  }else if(aggfunc=="medians"){
    aggfunc<-mediansByMask
  }else if(aggfunc=="custom"){
    aggfunc<-makeCustomAggregator(custom.aggfunc)
  }
  data<-as.data.frame(data)
  
  # check input
  if(check){
    datachecker(data=data,subjvar=subjvar,diffvars=diffvars,stratvars=stratvars,
                subscorevar=subscorevar,aggvar=aggvar,standardize=standardize)
  }
  
  # Arrange properly
  runorder<-do.call(order,data[c(subjvar,subscorevar,diffvars,stratvars)])
  arr.ds<- data[runorder,c(subjvar,subscorevar,diffvars,stratvars,aggvar)]
  origpps<-unique(arr.ds[[subjvar]])
  arr.ds[c(subjvar,subscorevar,diffvars,stratvars)]<-
    cols2ids(arr.ds[c(subjvar,subscorevar,diffvars,stratvars)])
  
  # Create indices
  arr.ds[[".diffidx"]]<-do.call(paste,arr.ds[,c(subjvar,subscorevar,diffvars),drop=FALSE])
  arr.ds[[".subscore"]]<-do.call(paste,arr.ds[,c(subjvar,subscorevar),drop=FALSE])
  subscores<-unique(arr.ds[[".subscore"]])
  pps<-unique(arr.ds[[subjvar]])
  
  diffidx<-unique(arr.ds[,c(subjvar,subscorevar,diffvars,".subscore",".diffidx"),drop=FALSE])
  subscorelist<-split(arr.ds,arr.ds[[".subscore"]])
  difflist<-split(arr.ds,arr.ds[[".diffidx"]])
  
  # Compute splithalves by subscore
  keys<-setNames(vector(mode="list",length=length(subscores)),subscores)
  if(verbose){
    cat("Generating splits...\n")
    pb = txtProgressBar(min = 0, max = length(subscores), initial = 0)
  }
  for(ss in subscores){
    if(verbose){ setTxtProgressBar(pb,which(ss==subscores)) }
    iterds<-arr.ds[arr.ds[[".subscore"]]==ss,c(diffvars,stratvars),drop=FALSE]
    iterrle<-rle(do.call(paste,args=iterds))
    grsizes<-iterrle$lengths
    if(length(grsizes)==0){grsizes<-nrow(iterds)}
    keys[[ss]]<-stratifiedItersplits(splits=splits, groupsizes=grsizes)
  }
  if(verbose){close(pb)}
  
  # get the other half by inverting the logical matrix
  antikeys<-lapply(keys,function(x) !x)
  
  # Get aggregates per condition
  keymeans<-antikeymeans<-matrix(ncol=splits,nrow=length(difflist))
  if(verbose){
    cat("Aggregating conditions...\n")
    pb = txtProgressBar(min = 0, max = length(difflist), initial = 0) 
  }
  rownames(keymeans)<-rownames(antikeymeans)<-names(difflist)
  for(i in seq_along(difflist)){
    if(verbose){ setTxtProgressBar(pb,i) }
    ssid<-difflist[[i]][[".subscore"]][1]
    diffid<-names(difflist)[i]
    keymeans[i,]<-
      aggfunc(x=difflist[[i]][[aggvar]],
              mask=keys[[ssid]][subscorelist[[ssid]][[".diffidx"]]==diffid,,drop=FALSE])
    antikeymeans[i,]<-
      aggfunc(x=difflist[[i]][[aggvar]],
              mask=antikeys[[ssid]][subscorelist[[ssid]][[".diffidx"]]==diffid,,drop=FALSE])
  }
  if(verbose){close(pb)}
  
  # Get scores
  if(length(diffvars)>0){
    diffidx[[".valence"]]<-Reduce(`*`,x=clamp.range(diffidx[diffvars]))
    keyscores<-antikeyscores<-matrix(ncol=splits,nrow=length(subscores))
    rownames(keyscores)<-rownames(antikeyscores)<-subscores
    for(i in seq_along(subscores)){
      keyscores[i,]<-
        colSums(keymeans[diffidx[[".subscore"]]==subscores[i],,drop=FALSE] * 
                  diffidx[[".valence"]][diffidx[[".subscore"]]==subscores[i]])
      
      antikeyscores[i,]<-
        colSums(antikeymeans[diffidx[[".subscore"]]==subscores[i],,drop=FALSE] * 
                  diffidx[[".valence"]][diffidx[[".subscore"]]==subscores[i]])
      
    }
  }else{
    keyscores<-keymeans
    antikeyscores<-antikeymeans
  }
  
  # Standardize if necessary
  if(standardize){
    keysds<-antikeysds<-matrix(ncol=splits,nrow=length(subscorelist))
    rownames(keysds)<-rownames(antikeysds)<-subscores
    if(verbose){
      cat("Computing standard deviations to standardize scores by...\n")
      pb = txtProgressBar(min = 0, max = length(subscorelist), initial = 0) 
    }
    for(ss in subscores){
      if(verbose){ setTxtProgressBar(pb,which(ss==subscores)) }
      keysds[ss,] <- sdsByMask(x=subscorelist[[ss]][[aggvar]], mask=keys[[ss]])
      antikeysds[ss,] <- sdsByMask(x=subscorelist[[ss]][[aggvar]], mask=antikeys[[ss]])
    }
    keyscores<-keyscores/keysds
    antikeyscores<-antikeyscores/antikeysds
    if(verbose){ close(pb) }
  }
  
  # Aggregate into single scores, if there were subscores
  if(length(subscorevar)>0){
    newkeyscores<-newantikeyscores<-matrix(ncol=splits,nrow=length(pps))
    rownames(newkeyscores)<-rownames(newantikeyscores)<-pps
    for(i in seq_along(pps)){
      ss<-unique(arr.ds[[".subscore"]][arr.ds[[subjvar]]==pps[i]])
      newkeyscores[i,]<-colMeans(keyscores[ss,])
      newantikeyscores[i,]<-colMeans(antikeyscores[ss,])
    }
    keyscores<-newkeyscores
    antikeyscores<-newantikeyscores
  }
  
  # Rename matrix rows to fit names in dataset
  rownames(keyscores)<-origpps[match(rownames(keyscores),pps)]
  rownames(antikeyscores)<-origpps[match(rownames(antikeyscores),pps)]
  
  # Get correlations
  cors<-corByColumns(keyscores,antikeyscores) |> SpearmanBrown()
  sampsize<-length(pps)
  
  # Form and return output
  out<-list(r=cormean(cors,sampsize),
            allcors=cors,
            nobs=sampsize)
  
  # Add individual split halves if requested
  if(include.scores){
    out$scores<-list(half1=keyscores,
                     half2=antikeyscores)
  }
  
  # Add class
  out<-structure(out,class="rapidsplit")
  
  # return
  return(out)
}



#' @param x \code{rapidsplit} object to print or plot.
#' @param ... Ignored.
#' @export
#' @rdname rapidsplit
print.rapidsplit<-function(x,...){
  coefstr<-paste0("Full-length reliability (Spearman-Brown coefficient):\n",
                  "SB (",mf(x$nobs),") = ",mf(x$r),
                  ", 95%CI [", mf(quantile(x$allcors,.025)), 
                  ", ", mf(quantile(x$allcors,.975)),"]",
                  ", p = ",mf(r2p(x$r,x$nobs),digits=3),"\n")
  cat(coefstr,sep="")
}

#' @param type Character argument indicating what should be plotted. 
#' By default, this plots the random split whose correlation is closest to the average.
#' However, this can also plot the random split with 
#' the \code{"minimum"} or \code{"maximum"} split-half correlation, or any \code{"random"} split. 
#' \code{"all"} splits can also be plotted together in one figure.
#' @param show.labels Should participant IDs be shown above their points in the scatterplot?
#' Defaults to \code{TRUE} and is ignored when \code{type} is \code{"all"}.
#'
#' @export
#' @rdname rapidsplit
plot.rapidsplit<-function(x,type=c("average","minimum","maximum","random","all"),
                          show.labels=TRUE,...){
  if(!("scores" %in% names(x))){
    message("Individual split-half scores were not included in rapidsplit object, so cannot be plotted")
    return()
  }
  
  type<-match.arg(type)
  if(type=="average"){
    title<-"Split-half Scatterplot for\nIteration with Average Reliability"
    idx<-which.min(abs(x$allcors-x$r))
  }else if(type=="minimum"){
    title<-"Split-half Scatterplot for\nIteration with the Lowest Reliability"
    idx<-which.min(x$allcors)
  }else if(type=="maximum"){
    title<-"Split-half Scatterplot for\nIteration with the Highest Reliability"
    idx<-which.max(x$allcors)
  }else if(type=="random"){
    title<-"Split-half Scatterplot for Random Iteration"
    idx<-sample(seq_along(x$allcors),1)
  }else if(type=="all"){
    title<-"Split-half Scatterplot for All Iterations"
    idx<-seq_along(x$allcors)
  }
  title <- paste0(title,"\n(r = ", mf(ifelse(length(idx)==1,x$allcors[idx],x$r)),")")
  
  h1vals<-x$scores$half1[,idx] |> as.vector()
  h2vals<-x$scores$half2[,idx] |> as.vector()
  plot(h1vals,h2vals,pch=20,main=title,
       xlab="Half 1 computed score",ylab="Half 2 computed score",
       col=rgb(0,0,0,ifelse(length(idx)==1,1,1/sqrt(length(idx)))))
  if(length(idx)==1 & show.labels){
    text(h1vals,h2vals,rownames(x$scores$half1),cex= 0.7, pos=3, offset=0.3)
  }
}



