
#include <Rcpp.h>
using namespace Rcpp;

//' colMedians
//' 
//' get column medians
//'
//' @param mat the matrix to retrieve column medians from
//' @export
// [[Rcpp::export]]
NumericVector colMedians(NumericMatrix mat){
  int tl = mat.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    NumericVector currcol = mat.column(i);
    out[i]=median(currcol);
  }
  return out;
}

//' @rdname colMedians
//' @param mat a matrix with values to aggregate
//' @param mask a logical matrix determining which data points to include and which not to
//' @export
// [[Rcpp::export]]
NumericVector colMedians_mask(NumericMatrix mat, LogicalMatrix mask){
  int tl = mat.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector currcol = mat.column(i);
    NumericVector maskedcol = currcol[currmask];
    out[i] = median(maskedcol);
  }
  return out;
}

//' @rdname colMedians
//' @param values Values to aggregate over in different mask configurations
//' @export
// [[Rcpp::export]]
NumericVector mediansByMask(NumericVector values, LogicalMatrix mask){
  int tl = mask.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector maskedcol = values[currmask];
    out[i] = median(maskedcol);
  }
  return out;
}

//' @rdname colMedians
//' @export
// [[Rcpp::export]]
NumericVector colMeans_mask(NumericMatrix mat, LogicalMatrix mask){
  int tl = mat.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector currcol = mat.column(i);
    NumericVector maskedcol = currcol[currmask];
    out[i] = mean(maskedcol);
  }
  return out;
}

//' @rdname colMedians
//' @param values Values to aggregate over in different mask configurations
//' @export
// [[Rcpp::export]]
NumericVector meansByMask(NumericVector values, LogicalMatrix mask){
  int tl = mask.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector maskedcol = values[currmask];
    out[i] = mean(maskedcol);
  }
  return out;
}

//' colSds
//' 
//' get column SDs
//'
//' @param mat the matrix to retrieve column SDs from
//' @export
// [[Rcpp::export]]
NumericVector colSds(NumericMatrix mat){
  int tl = mat.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    NumericVector currcol = mat.column(i);
    out[i]=sd(currcol);
  }
  return out;
}

//' @rdname colSds
//' @param mask a logical matrix determining which data points to include and which not to
//' @export
// [[Rcpp::export]]
NumericVector colSds_mask(NumericMatrix mat, LogicalMatrix mask){
  int tl = mat.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector currcol = mat.column(i);
    NumericVector maskedcol = currcol[currmask];
    out[i]=sd(maskedcol);
  }
  return out;
}

//' @rdname colSds
//' @param values Values to aggregate over in different mask configurations
//' @export
// [[Rcpp::export]]
NumericVector sdsByMask(NumericVector values, LogicalMatrix mask){
  int tl = mask.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector maskedcol = values[currmask];
    out[i]=sd(maskedcol);
  }
  return out;
}

//' Exclude SD-based outliers
//' 
//' Update a mask matrix based on outlyingness
//' 
//' @param rtvec Reaction time vector
//' @param mask a logical matrix determining which data points to include and which not to
//' @param sdlim Standard deviation limit to apply; values beyond are classified as outliers and masked
//' @returns An updated mask
//' @export
// [[Rcpp::export]]
LogicalMatrix ExcludeSDOutliers(NumericVector rtvec, LogicalMatrix mask, double sdlim = 3){
  
  int outcols = mask.ncol();
  
  LogicalMatrix newmask = mask;
  
  for(int i=0; i<outcols; i++){
    LogicalVector currmask = mask.column(i);
    NumericVector maskedcol = rtvec[currmask];
    
    double colmean = mean(maskedcol);
    double colsd = sd(maskedcol);
    
    LogicalVector newmaskedmask = 
      (maskedcol > (colmean - sdlim * colsd)) & 
      (maskedcol < (colmean + sdlim * colsd));
    
    currmask[currmask] = newmaskedmask;
    newmask(_,i) = currmask;
  }
  
  return newmask;
}

//' @rdname ExcludeSDOutliers
//' @param mat Matrix in which to mark SD-based outleirs by column (with FALSE)
//' @export
// [[Rcpp::export]]
LogicalMatrix ExcludeSDOutliers_nomask(NumericMatrix mat, double sdlim = 3){
  
  int outrows = mat.nrow();
  int outcols = mat.ncol();
  
  LogicalMatrix newmask (outrows,outcols);
  
  for(int i=0; i<outcols; i++){
    NumericVector currcol = mat.column(i);
    
    double colmean = mean(currcol);
    double colsd = sd(currcol);
    
    LogicalVector currnewmask = 
      (currcol > (colmean - sdlim * colsd)) & 
      (currcol < (colmean + sdlim * colsd));
    
    newmask(_,i) = currnewmask;
  }
  
  return newmask;
}

IntegerMatrix itersplits(int x, IntegerVector y, bool replace=false){
  int ysize = y.size();
  IntegerMatrix itermat (ysize,x);
  for(int i = 0; i < x; i++) {
    itermat(_,i) = Rcpp::sample(y, ysize, replace);
  }
  return itermat;
}

//' applyItersplits
//' 
//' generate splits for splithalf
//'
//' @param iters number of iterations
//' @param splits list of vectors of row numbers
//' @param replace Sample without (default) or with replacement
//' @export
// [[Rcpp::export]]
List applyItersplits(int iters, List splits, bool replace=false){
  int tl = splits.size();
  List out (tl);
  
  for(int i=0; i<tl; i++){
    out[i]=itersplits(iters,splits[i],replace);
  }
  return out;
}

//' stratified_itersplits
//' 
//' generate stratified splits for a single participant
//' 
//' This first equally splits what can be equally split within groups.
//' Then it randomly splits all the leftovers.
//'
//' @param itercount number of iterations
//' @param groupsizes vector of number of RTs per group to stratify
//' @returns A matrix with zeroes and ones
//' @export
// [[Rcpp::export]]
IntegerMatrix stratified_itersplits(int itercount, IntegerVector groupsizes){
  
  // Get info on arguments and create output matrix
  int sampsize = sum(groupsizes);
  int groupcount = groupsizes.length();
  IntegerMatrix itermat (sampsize,itercount);
  
  // count uneven subgroups and prepare an empty vector 
  // to hold half assignments for their leftovers
  int n_uneven = 0;
  for(int i = 0; i < groupcount; i++){
    n_uneven += groupsizes[i] % 2;
  }
  IntegerVector uneven_sampvec(n_uneven);
  
  // Create simple vectors for use in assignment to halves
  IntegerVector zeroone = {0,1};
  IntegerVector onezero = {1,0};
  IntegerVector sampvec_zeroone = rep_len(zeroone,n_uneven);
  IntegerVector sampvec_onezero = rep_len(onezero,n_uneven);
  
  // loop over iterations (columns in output matrix)
  for(int i = 0; i < itercount; i++) {
    
    // Generate an assignments vector for all the leftovers 
    // after all evenly splittable values have been split
    bool sampvec_order = sample(onezero,1,false);
    if(sampvec_order){
      uneven_sampvec = sample(sampvec_onezero,n_uneven,false);
    }else{
      uneven_sampvec = sample(sampvec_zeroone,n_uneven,false);
    }
    
    // Assign all to halves
    int min_index = 0;
    int uneven_index = 0;
    for(int j = 0; j < groupcount; j++){
      IntegerVector sampvec = rep_len(zeroone, groupsizes[j]);
      if(groupsizes[j] % 2 == 1){
        sampvec[0] = uneven_sampvec[uneven_index];
        uneven_index++;
      }
      sampvec = Rcpp::sample(sampvec, groupsizes[j], false);
      
      for(int k = 0; k < groupsizes[j]; k++){
        itermat(min_index+k,i)=sampvec[k];
      }
      min_index += groupsizes[j];
    }
  }
  
  return itermat;
}


//' Correlate each column of 1 matrix with the same column in another matrix
//' 
//' Correlate each column of 1 matrix with the same column in another matrix
//'
//' @param mat1,mat2 Matrices whose values to correlate by column
//' @returns A numeric vector of correlations per column
//' @export
// [[Rcpp::export]]
NumericVector corByColumns(NumericMatrix mat1, NumericMatrix mat2){
  int tl = mat1.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    NumericVector currmat1 = mat1.column(i);
    NumericVector currmat2 = mat2.column(i);
    
    LogicalVector incl = (!is_na(currmat1)) & (!is_na(currmat2));
    currmat1 = currmat1[incl];
    currmat2 = currmat2[incl];
    double currn = sum(incl);
    
    out[i] = (1/(currn-1)) *
      sum( ((currmat1-mean(currmat1))/sd(currmat1)) * 
      ((currmat2-mean(currmat2))/sd(currmat2)) );
  }
  return out;
}

//' @rdname corByColumns
//' @param mask Logical matrix marking which data points to include
//' @export
// [[Rcpp::export]]
NumericVector corByColumns_mask(NumericMatrix mat1, NumericMatrix mat2, LogicalMatrix mask){
  int tl = mask.ncol();
  NumericVector out (tl);
  
  for(int i=0; i<tl; i++){
    NumericVector currmat1 = mat1.column(i);
    NumericVector currmat2 = mat2.column(i);
    
    LogicalVector incl = mask.column(i);
    currmat1 = currmat1[incl];
    currmat2 = currmat2[incl];
    double currn = sum(incl);
    
    out[i] = (1/(currn-1)) *
      sum( ((currmat1-mean(currmat1))/sd(currmat1)) * 
      ((currmat2-mean(currmat2))/sd(currmat2)) );
  }
  return out;
}

