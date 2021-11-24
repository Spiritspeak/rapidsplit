

#' rapidsplit
#' ultra fast split-half
#'
#' @param ds dataset
#' @param subjvar subject var name
#' @param pullvar movement direction var name
#' @param targetvar stim type var name
#' @param rtvar rt varname
#' @param iters n iterations
#' @param agg means or medians
#' @param standardize divide by sd or not
#'
#' @return
#' @export
#'
#' @examples
#' print("no example")
rapidsplit<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,agg=c("means","medians"),standardize=F){
  if(agg=="means"){
    aggfunc<-colMeans
  }else if(agg=="medians"){
    aggfunc<-function(x){apply(x,2,median.default)}
  }
  
  presplit<-split(seq_len(nrow(ds)),ds[c(subjvar,pullvar,targetvar)])
  
  origkeys<-applyItersplits(iters,presplit)
  keys<-lapply(origkeys,function(x){ x[1:floor(nrow(x)/2),] } )
  antikeys<-lapply(origkeys,function(x)x[(floor(nrow(x)/2)+1):nrow(x),])
  rm(origkeys)
  
  keymeans<-t(sapply(keys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  antikeymeans<-t(sapply(antikeys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  
  sk<-gsub("^[0-9]*.","",names(presplit))
  
  keyscores<-(keymeans[sk=="0.0",] + keymeans[sk=="1.1",]) - (keymeans[sk=="0.1",] + keymeans[sk=="1.0",])
  antikeyscores<-(antikeymeans[sk=="0.0",] + antikeymeans[sk=="1.1",]) - (antikeymeans[sk=="0.1",] + antikeymeans[sk=="1.0",])
  
  if(standardize){
    sj<-gsub("\\..*$","",names(presplit))
    
    subjkeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,keys[x==sj])),unique(sj))
    keysds<-t(sapply(subjkeys,function(x) apply(matrix(ds[[rtvar]][x],ncol=iters),2,sd)))
    keyscores<-keyscores/keysds
    
    subjantikeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,antikeys[x==sj])),unique(sj))
    antikeysds<-t(sapply(subjantikeys,function(x) apply(matrix(ds[[rtvar]][x],ncol=iters),2,sd)))
    antikeyscores<-antikeyscores/antikeysds
  }
  
  cors<-sapply(1:iters,function(x) cor(keyscores[,x],antikeyscores[,x]))
  
  list(r=cormean(SpearmanBrown(cors),rep(length(unique(ds[[subjvar]])),length(cors))),
       allcors=cors)
}