#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix cpp_gensnp(const int n, const int p, const float maf=0.5){

  float step1 = pow((1-maf),2);
  float step2 = pow((1-maf),2) + 2*maf*(1-maf);
  NumericMatrix unif(n, p);
  
  for(int i = 0; i< n; i++){
    for(int j = 0; j < p; j++){
       float a = (float) rand()/RAND_MAX;
       if(a > step2)
         unif(i,j) = 2;
       else if(a <= step2 & a > step1)
         unif(i,j) = 1;
       else
         unif(i,j) = 0;
    }
  }
  return(unif);
}










