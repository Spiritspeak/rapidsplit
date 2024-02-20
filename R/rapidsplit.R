
cols2ids<-function(df){
  for(col in seq_along(df)){
    df[,col]<-as.numeric(as.factor(df[,col]))
  }
  df
}

clamp.range<-function(x){
  x[x==min(x)]<- -1
  x[x==max(x)]<- 1
  x
}

#' rapidsplit
#' @description A very fast algorithm for permutated split-half reliability
#'
#' @param ds Dataset, a \code{data.frame}
#' @param subjvar Subject ID variable name, a \code{character}
#' @param diffvars Variables that determine which conditions 
#' need to be subtracted from each other, a \code{character}
#' @param stratvars Additional variables that the splits should 
#' be stratified by, if possible; a \code{character}
#' @param rtvar Reaction time variable name, a \code{character}
#' @param iters Number of split-halves to average, an \code{integer}
#' @param agg The function by which to aggregate the RTs; can be \code{"means"} or \code{"medians"}
#' @param standardize Whether to divide by scores by the subject's SD; a \code{logical}
#'
#' @return A \code{list} containing the averaged reliability as well as 
#' a vector with the reliability of each iteration
#' @export
#'
#' @examples
#' print("no example")
rapidsplit<-function(ds,subjvar,diffvars=NULL,stratvars=NULL,rtvar,iters,
                     agg=c("means","medians"),standardize=F){
  agg<-match.arg(agg)
  if(agg=="means"){
    aggfunc<-meansByMask
  }else if(agg=="medians"){
    aggfunc<-mediansByMask
  }
  
  # Arrange properly
  arr.ds<- ds[do.call(order,ds[c(subjvar,diffvars,stratvars)]),
              c(subjvar,diffvars,rtvar,stratvars)]
  arr.ds[,c(subjvar,diffvars,stratvars)] <-
    cols2ids(arr.ds[,c(subjvar,diffvars,stratvars),drop=F])
  arr.ds[[".diffidx"]]<-interaction(arr.ds[c(subjvar,diffvars)],lex.order=T)
  
  pps<-unique(arr.ds[,subjvar])
  diffidx<-unique(arr.ds[c(subjvar,diffvars,".diffidx")])
  subjlist<-split(arr.ds,arr.ds[[subjvar]])
  difflist<-split(arr.ds,arr.ds[[".diffidx"]])
  
  # Stratify by subject
  origkeys<-setNames(vector(mode="list",length=length(pps)),pps)
  for(sj in pps){
    iterds<-arr.ds[arr.ds[[subjvar]]==sj,c(diffvars,stratvars),drop=F]
    iterrle<-rle(do.call(paste,args=cols2ids(iterds)))
    grsizes<-iterrle$lengths
    if(length(grsizes)==0){grsizes<-nrow(iterds)}
    origkeys[[sj]]<-stratified_itersplits(itercount=iters, groupsizes=grsizes)
  }
  
  #take two halves of the data separately
  keys<-lapply(origkeys,function(x) x==1)
  antikeys<-lapply(origkeys,function(x) x==0)
  
  #free memory (these matrices get huge)
  rm(origkeys)
  
  #get aggregates per condition
  keymeans<-antikeymeans<-matrix(ncol=iters,nrow=length(difflist))
  for(i in seq_along(difflist)){
    ppid<-difflist[[i]][[subjvar]][1]
    diffid<-names(difflist)[i]
    keymeans[i,]<-
      aggfunc(values=difflist[[i]][[rtvar]],
              mask=keys[[ppid]][subjlist[[ppid]][[".diffidx"]]==diffid ,,drop=F])
    antikeymeans[i,]<-
      aggfunc(values=difflist[[i]][[rtvar]],
              mask=antikeys[[ppid]][subjlist[[ppid]][[".diffidx"]]==diffid ,,drop=F])
  }
  
  # Get scores
  if(length(diffvars)>0){
    diffidx[[".valence"]]<-Reduce(`*`,x=clamp.range(diffidx[diffvars]))
    keyscores<-t(sapply(pps,function(pp){
      colSums(keymeans[diffidx[[subjvar]]==pp,] * diffidx[[".valence"]][diffidx[[subjvar]]==pp])
    }))
    antikeyscores<-t(sapply(pps,function(pp){
      colSums(antikeymeans[diffidx[[subjvar]]==pp,] * diffidx[[".valence"]][diffidx[[subjvar]]==pp])
    }))
  }else{
    keyscores<-keymeans
    antikeyscores<-antikeymeans
  }
    
  # Standardize if necessary
  if(standardize){
    keysds<-antikeysds<-matrix(ncol=iters,nrow=length(subjlist))
    for(i in seq_along(subjlist)){
      keysds[i,]<-sdsByMask(values=subjlist[[i]][[rtvar]],
                            mask=keys[[i]])
      antikeysds[i,]<-sdsByMask(values=subjlist[[i]][[rtvar]],
                                mask=antikeys[[i]])
    }
    keyscores<-keyscores/keysds
    antikeyscores<-antikeyscores/antikeysds
  }
  
  # Get correlations
  cors<-corByColumns(keyscores,antikeyscores)
  
  # Form and return output
  out<-list(r=cormean(SpearmanBrown(cors),rep(length(unique(ds[[subjvar]])),length(cors))),
            allcors=cors)
  return(out)
}

