

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
    aggfunc<-colMedians
  }
  
  #Stratify the split
  presplit<-split(seq_len(nrow(ds)),ds[c(subjvar,pullvar,targetvar)])
  
  #obtain the splitting keys in a list of matrices, where each column is an iteration
  origkeys<-applyItersplits(iters,presplit)
  
  #take the upper and lower halves of each split matrix, thereby sampling two halves of the data
  keys<-lapply(origkeys,function(x){ x[1:floor(nrow(x)/2),] } )
  antikeys<-lapply(origkeys,function(x)x[(floor(nrow(x)/2)+1):nrow(x),])
  
  #free memory (these matrices get huge)
  rm(origkeys)
  
  keymeans<-t(sapply(keys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  antikeymeans<-t(sapply(antikeys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  
  sk<-gsub("^[0-9]*.","",names(presplit))
  
  keyscores<-(keymeans[sk=="0.0",] + keymeans[sk=="1.1",]) - (keymeans[sk=="0.1",] + keymeans[sk=="1.0",])
  antikeyscores<-(antikeymeans[sk=="0.0",] + antikeymeans[sk=="1.1",]) - (antikeymeans[sk=="0.1",] + antikeymeans[sk=="1.0",])
  
  if(standardize){
    sj<-gsub("\\..*$","",names(presplit))
    
    subjkeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,keys[x==sj])),unique(sj))
    keysds<-t(sapply(subjkeys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    keyscores<-keyscores/keysds
    
    subjantikeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,antikeys[x==sj])),unique(sj))
    antikeysds<-t(sapply(subjantikeys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    antikeyscores<-antikeyscores/antikeysds
  }
  
  cors<-sapply(1:iters,function(x) cor(keyscores[,x],antikeyscores[,x]))
  
  list(r=cormean(SpearmanBrown(cors),rep(length(unique(ds[[subjvar]])),length(cors))),
       allcors=cors)
}










stratRapidsplit2<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,agg=c("means","medians"),standardize=F,stratvars=NULL){
  if(agg=="means"){
    aggfunc<-colMeans
  }else if(agg=="medians"){
    aggfunc<-colMedians
  }
  
  #Stratify the split
  presplit<-split(seq_len(nrow(ds)),ds[c(subjvar,pullvar,targetvar,stratvars)],drop=T)
  
  #obtain the splitting keys in a list of matrices, where each column is an iteration
  origkeys<-applyItersplits(iters,presplit)
  
  #take the upper and lower halves of each split matrix, thereby sampling two halves of the data
  keys<-lapply(origkeys,function(x){ x[1:floor(nrow(x)/2),] } )
  antikeys<-lapply(origkeys,function(x)x[(floor(nrow(x)/2)+1):nrow(x),])
  
  #free memory (these matrices get huge)
  rm(origkeys)
  
  #remerge across stratification variables, if there were any
  if(length(stratvars>0)){
    sk<-gsub("\\.[^.]+$","",names(presplit))
    keys<-lapply(tapply(keys,sk,list),function(x)do.call(rbind,x))
    antikeys<-lapply(tapply(antikeys,sk,list),function(x)do.call(rbind,x))
  }
  
  #get aggregates per condition
  keymeans<-t(sapply(keys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  antikeymeans<-t(sapply(antikeys,function(x) aggfunc(matrix(ds[[rtvar]][x],ncol=iters))))
  
  sk<-gsub("^[0-9]*.","",names(keys))
  
  keyscores<-(keymeans[sk=="0.0",] + keymeans[sk=="1.1",]) - (keymeans[sk=="0.1",] + keymeans[sk=="1.0",])
  antikeyscores<-(antikeymeans[sk=="0.0",] + antikeymeans[sk=="1.1",]) - (antikeymeans[sk=="0.1",] + antikeymeans[sk=="1.0",])
  
  if(standardize){
    sj<-gsub("\\..*$","",names(presplit))
    
    subjkeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,keys[x==sj])),unique(sj))
    keysds<-t(sapply(subjkeys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    keyscores<-keyscores/keysds
    
    subjantikeys<-setNames(lapply(unique(sj),function(x)do.call(rbind,antikeys[x==sj])),unique(sj))
    antikeysds<-t(sapply(subjantikeys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    antikeyscores<-antikeyscores/antikeysds
  }
  
  cors<-sapply(1:iters,function(x) cor(keyscores[,x],antikeyscores[,x]))
  
  list(r=cormean(SpearmanBrown(cors),rep(length(unique(ds[[subjvar]])),length(cors))),
       allcors=cors)
}


stratRapidsplit<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,standardize=F,stratvars=NULL){
  
  #Stratify the split
  presplit<-split(seq_len(nrow(ds)),ds[c(subjvar,stratvars)],drop=T)
  
  #obtain the splitting keys in a list of matrices, where each column is an iteration
  origkeys<-applyItersplits(iters,presplit)
  
  #take the upper and lower halves of each split matrix, thereby sampling two halves of the data
  keys<-lapply(origkeys,function(x){ x[1:floor(nrow(x)/2),] } )
  antikeys<-lapply(origkeys,function(x)x[(floor(nrow(x)/2)+1):nrow(x),])
  
  #free memory (these matrices get huge)
  rm(origkeys)
  
  #remerge across stratifications
  if(length(stratvars>0)){
    sk<-gsub("\\..*$","",names(presplit))
    keys<-lapply(tapply(keys,sk,list),function(x)do.call(rbind,x))
    antikeys<-lapply(tapply(antikeys,sk,list),function(x)do.call(rbind,x))
  }
  
  
  keypulltar<-t(sapply(keys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==1 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==1) 
  }))
  keypushtar<-t(sapply(keys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==0 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==1) 
  }))
  keypullcon<-t(sapply(keys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==1 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==0) 
  }))
  keypushcon<-t(sapply(keys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==0 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==0) 
  }))
  keyscores<-((keypushtar-keypulltar)-(keypushcon-keypullcon))
  
  
  antikeypulltar<-t(sapply(antikeys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==1 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==1) 
  }))
  antikeypushtar<-t(sapply(antikeys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==0 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==1) 
  }))
  antikeypullcon<-t(sapply(antikeys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==1 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==0) 
  }))
  antikeypushcon<-t(sapply(antikeys,function(x){ 
    colMeans_mask(mat = matrix(ds[[rtvar]][x],ncol=iters),
                  mask= matrix(ds[[pullvar]][x],ncol=iters)==0 & 
                    matrix(ds[[targetvar]][x],ncol=iters)==0) 
  }))
  antikeyscores<-((antikeypushtar-antikeypulltar)-(antikeypushcon-antikeypullcon))
  
  # divide by SD if necessary
  if(standardize){
    keysds<-t(sapply(keys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    keyscores<-keyscores/keysds
    
    antikeysds<-t(sapply(antikeys,function(x) colSds(matrix(ds[[rtvar]][x],ncol=iters))))
    antikeyscores<-antikeyscores/antikeysds
  }
  
  cors<-sapply(1:iters,function(x) cor(keyscores[,x],antikeyscores[,x]))
  
  list(r=cormean(SpearmanBrown(cors),rep(length(unique(ds[[subjvar]])),length(cors))),
       allcors=cors)
}

