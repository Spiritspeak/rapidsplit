
#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix itersplits(int x, NumericVector y){
  int ysize = y.size();
  NumericMatrix itermat (ysize,x);
  for(int i = 0; i < x; i++) {
    itermat(_,i) = Rcpp::sample(y, ysize);
  }
  return itermat;
}

//' applyItersplits
//' 
//' generate splits for splithalf
//'
//' @param iters number of iterations
//' @param splits list of vectors of row numbers
//' @export
// [[Rcpp::export]]
List applyItersplits(int iters, List splits){
  int tl = splits.size();
  List out (tl);
  
  for(int i=0; i<tl; i++){
    out[i]=itersplits(iters,splits[i]);
  }
  return out;
}