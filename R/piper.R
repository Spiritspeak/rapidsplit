#data(foodAAT)

generateSplits<-function(data,subsetvars,stratvars=NULL,splits,verbose=T){
  
  # Arrange properly
  runorder<-do.call(order,data[,c(subsetvars,stratvars),drop=FALSE])
  arr.ds<- cols2ids(data[runorder,c(subsetvars,stratvars),drop=FALSE])
  
  subsetvec<-do.call(paste,arr.ds[,subsetvars,drop=FALSE])
  subsets<-unique(subsetvec)
  
  # Stratify by split stratum
  origkeys<-setNames(vector(mode="list",length=length(subsets)),subsets)
  if(verbose){ pb = txtProgressBar(min = 0, max = length(subsets), initial = 0) }
  for(ss in subsets){
    if(verbose){ setTxtProgressBar(pb,which(ss==subsets)) }
    iterds<-arr.ds[subsetvec==ss,c(stratvars),drop=FALSE]
    iterrle<-rle(do.call(paste,args=iterds))
    grsizes<-iterrle$lengths
    if(length(grsizes)==0){grsizes<-nrow(iterds)}
    origkeys[[ss]]<-stratifiedItersplits(splits=splits, groupsizes=grsizes)
  }
  origkeys<-do.call(rbind,origkeys)
  if(verbose){ close(pb) }
  
  # return
  backorder<-order(runorder)
  list(indices=arr.ds[backorder,c(subsetvars,stratvars)],
       keys=origkeys[backorder,])
}

# test<-generateSplits(data=foodAAT,subjvar="subjectid",
#                      stratvars=c("is_pull","is_target","stimid"),splits=5000)

applyAggregator<-function(data,subsetvars,aggvar,aggfunc,mask,verbose=T){
  subsetvec<-do.call(paste,cols2ids(data[,subsetvars,drop=FALSE]))
  runorder<-order(subsetvec)
  subsetvec<-subsetvec[runorder]
  data<-data[runorder,]
  mask<-mask[runorder,]
  
  subsets<-unique(subsetvec)
  aggs<-matrix(ncol=ncol(mask),nrow=length(subsets))
  if(verbose){ pb = txtProgressBar(min = 0, max = length(subsets), initial = 0) }
  for(i in seq_along(subsets)){
    if(verbose){ setTxtProgressBar(pb,i) }
    aggs[i,]<-
      aggfunc(x=data[[aggvar]][subsetvec==subsets[i]],
              mask=mask[subsetvec==subsets[i],,drop=F])
  }
  if(verbose){ close(pb) }
  list(indices=unique(data[,subsetvars]),aggs=aggs)
}

applyAggregator2<-function(data,subsetvars,aggvar,aggfunc,masks,verbose=T){
  outputs<-list()
  for(i in seq_along(masks)){
    outputs[[i]]<-applyAggregator(data,subsetvars,aggvar,aggfunc,masks[[i]],verbose)
  }
  outputs
}

# stuff<-applyAggregator(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                        aggfunc=meansByMask,mask=test$keys)

reduceRows<-function(indices,aggs,groupvars,collapsevar,
                     action=c("subtract","add","divide","multiply"),
                     highlow=T){
  action<-match.arg(action)
  
  subsetvec<-do.call(paste,indices[,c(groupvars,collapsevar),drop=FALSE])
  if(is.unsorted(subsetvec)){
    runorder<-order(subsetvec)
    indices<-indices[runorder,]
    aggs<-aggs[runorder,]
  }
  
  groupvec<-do.call(paste,indices[,c(groupvars),drop=FALSE])
  
  # add
  if(action=="add" | action=="multiply"){
    subsets<-unique(groupvec)
    newagg<-matrix(ncol=ncol(aggs),nrow=length(subsets))
    agger<-ifelse(action=="add",colSums,colProds)
    for(i in seq_along(subsets)){
      newagg[i,]<-agger(aggs[groupvec==subsets[i],])
    }
    newindices<-unique(indices[,c(groupvars),drop=FALSE])
  }else{
    # Checks
    condcounts<-table(groupvec)
    if(!all(sapply(unique(groupvec),
                   function(x){length(unique(indices[groupvec==x,collapsevar]))}) == 2)){
      stop("More or less than 2 levels of collapse variable detected.")
    }
    
    if(any(condcounts!=2)){
      stop("In ",sum(condcounts!=2)," case(s) there are less or more than 2 instances",
           " of a permutation of grouping variables.",
           " Cannot unambiguously match the collapse variable.")
    }
    if(highlow){
      mat1<-aggs[indices[[collapsevar]]==max(indices[[collapsevar]]),]
      mat2<-aggs[indices[[collapsevar]]==min(indices[[collapsevar]]),]
    }else{
      mat1<-aggs[indices[[collapsevar]]==min(indices[[collapsevar]]),]
      mat2<-aggs[indices[[collapsevar]]==max(indices[[collapsevar]]),]
    }
    if(action=="subtract"){ newagg<-mat1-mat2 }else
    if(action=="divide"){ newagg<-mat1/mat2 }
    newindices<-indices[indices[[collapsevar]]==max(indices[[collapsevar]]),
                        colnames(indices)!=collapsevar]
  }
  list(indices=newindices,aggs=newagg)
}

# lessstuff<-
#   reduceRows(indices=stuff$indices,aggs=stuff$aggs,groupvars=c("subjectid","is_target"),
#            collapsevar="is_pull",highlow=F)
# lesserstuff<-
#   reduceRows(indices=lessstuff$indices,aggs=lessstuff$aggs,groupvars=c("subjectid"),
#              collapsevar="is_target")
# 
# 

# library(tictoc)
# library(rapidsplit)
# data(foodAAT)
# tic()
# splits<-generateSplits(data=foodAAT,subsetvars="subjectid",
#                        stratvars=c("is_pull","is_target","stimid"),splits=5000)
# 
# aggs<-applyAggregator(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                       aggfunc=meansByMask,mask=splits$keys)
# lessstuff<-
#   reduceRows(indices=aggs$indices,aggs=aggs$aggs,groupvars=c("subjectid","is_target"),
#              collapsevar="is_pull",highlow=F)
# lesserstuff<-
#   reduceRows(indices=lessstuff$indices,aggs=lessstuff$aggs,groupvars=c("subjectid"),
#              collapsevar="is_target")
# 
# aggs2<-applyAggregator(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                       aggfunc=meansByMask,mask=!splits$keys)
# lessstuff2<-
#   reduceRows(indices=aggs2$indices,aggs=aggs2$aggs,groupvars=c("subjectid","is_target"),
#              collapsevar="is_pull",highlow=F)
# lesserstuff2<-
#   reduceRows(indices=lessstuff2$indices,aggs=lessstuff2$aggs,groupvars=c("subjectid"),
#              collapsevar="is_target")
# 
# allcors<-corByColumns(lesserstuff$aggs,lesserstuff2$aggs)
# SpearmanBrown(allcors) |> cormean(36)
# toc()
# 
# tic()
# rapidsplit(data=foodAAT,subjvar="subjectid",diffvars=c("is_pull","is_target"),
#            stratvars=c("stimid"),aggvar="RT",splits=5000,standardize=F)
# toc()
# 
# 
# lesserstuff2<-
#   reduceRows(indices=lessstuff2$indices,aggs=lessstuff2$aggs,groupvars=c("subjectid"),
#              collapsevar="is_target",action="add")
# 
# # Allow multiple aggs as input (list)
# # put data input vars first
# 
# 
# h<-applyAggregator(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                    aggfunc=iqragger,mask=splits$keys)
# h2<-applyAggregator(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                    aggfunc=sdsByMask,mask=splits$keys)
# plot(h$aggs[,1],h2$aggs[,1])
# 
# splits2<-generateSplits2(data=foodAAT,splitvars=c("subjectid","is_pull","is_target"),
#                          stratvars=c("stimid"),splits=5000)
# 
# 
# h<-as.numeric(splits$keys)
# 
# all.equal(h==1,as.logical(h))
# 
# microbenchmark::microbenchmark(h==1,as.logical(h))
# 
# 
# hh<-applyAggregator2(foodAAT,subsetvars=c("subjectid","is_pull","is_target"),aggvar="RT",
#                      aggfunc=meansByMask,mask=list(splits$keys,!splits$keys))
# 
# rapidsplit(data=foodAAT,subjvar="subjectid",diffvars=c("is_pull"),stratvars=c("stimid"),
#             aggvar="RT",splits=1000)




