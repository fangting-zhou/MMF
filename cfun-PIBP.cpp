// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
double updater(mat B, double a = 1, double b = 1) {
  a = a + sum(vectorise(B));
  b = b + sum(vectorise(1 - B));
  
  return(R::rbeta(a, b));
}

// [[Rcpp::export]]
mat updateB(mat Z, mat A, mat B, mat W, vec z, double r) {
  int i, j, k; mat BNEW = B;
  
  for(j = 0; j < B.n_rows; j ++) {
    for(k = 0; k < B.n_cols; k ++) {
      BNEW(j, k) = 1 - B(j, k);
      
      double postold = log(r) * B(j, k) + log(1 - r) * (1 - B(j, k));
      double postnew = log(r) * BNEW(j, k) + log(1 - r) * (1 - BNEW(j, k));
      
      for(i = 0; i < A.n_rows; i ++) {
        double oddsold = sum(A.row(i) % W.row(i) % B.row(j)) + z(i);
        double oddsnew = sum(A.row(i) % W.row(i) % BNEW.row(j)) + z(i);
        
        postold += oddsold * Z(i, j) - log(1 + exp(oddsold));
        postnew += oddsnew * Z(i, j) - log(1 + exp(oddsnew));
      }
      
      if((postold - postnew) <= log(1 / R::runif(0, 1) - 1)) B(j, k) = BNEW(j, k); else BNEW(j, k) = B(j, k);
    }
  }
  
  return(B);
}

// [[Rcpp::export]]
mat updateW(mat Z, mat A, mat B, mat W, vec z, double sw = 0.1, double s = 1) {
  int i, k; mat WNEW = W;
  
  for(i = 0; i < W.n_rows; i ++) {
    for(k = 0; k < W.n_cols; k ++) {
      do {
        WNEW(i, k) = R::rnorm(W(i, k), s);
      } while (WNEW(i, k) <= 0);
      
      vec oddsnum = B * (A.row(i) % WNEW.row(i)).t() + z(i);
      vec oddsden = B * (A.row(i) % W.row(i)).t() + z(i);
      
      double num = dot(oddsnum, Z.row(i)) - sum(log(1 + exp(oddsnum))) - sw * WNEW(i, k);
      double den = dot(oddsden, Z.row(i)) - sum(log(1 + exp(oddsden))) - sw * W(i, k);
      
      if(log(R::runif(0, 1)) <= (num - den)) W(i, k) = WNEW(i, k); else WNEW(i, k) = W(i, k);
    }
  }
  
  return(W);
}

// [[Rcpp::export]]
vec updatez(mat Z, mat A, mat B, mat W, vec z, double sz = 100, double s = 1) {
  int i;
  
  for(i = 0; i < z.n_elem; i ++) {
    double znew = R::rnorm(z(i), s);
    
    vec oddsnum = B * (A.row(i) % W.row(i)).t() + znew;
    vec oddsden = B * (A.row(i) % W.row(i)).t() + z(i);
      
    double num = - pow(znew, 2) / pow(sz, 2) / 2 + dot(oddsnum, Z.row(i)) - sum(log(1 + exp(oddsnum)));
    double den = - pow(z(i), 2) / pow(sz, 2) / 2 + dot(oddsden, Z.row(i)) - sum(log(1 + exp(oddsden)));
    
    if(R::runif(0, 1) <= (num - den)) z(i) = znew;
  }
  
  return(z);
}

// [[Rcpp::export]]
mat updateZ(mat X, mat Z, mat A, mat B, mat W, vec z, vec a, vec b) {
  int i, j;
  
  for(i = 0; i < Z.n_rows; i ++) {
    for(j = 0; j < Z.n_cols; j ++) {
      double prob = sum(A.row(i) % W.row(i) % B.row(j)) + z(i);
      double plus = sum(b + (a - b) % Z.col(j)) - b(i) - (a(i) - b(i)) * Z(i, j);
      
      double postone = lgamma(plus + a(i)) - lgamma(sum(X.col(j)) + plus + a(i)) +
        lgamma(X(i, j) + a(i)) - lgamma(a(i)) + prob;
      double postzero = lgamma(plus + b(i)) - lgamma(sum(X.col(j)) + plus + b(i)) +
        lgamma(X(i, j) + b(i)) - lgamma(b(i));
      
      Z(i, j) = ((postzero - postone) <= log(1 / R::runif(0, 1) - 1));
    }
  }
  
  return(Z);
}

// [[Rcpp::export]]
vec updatea(mat X, mat Z, vec a, vec b, double alpha = 1, double beta = 0.1, double s = 1) {
  int i, j; vec anew = a;
  
  for(i = 0; i < a.n_elem; i ++) {
    do {
      anew(i) = R::rnorm(a(i), s);
    } while (anew(i) <= b(i));
    
    double num = (alpha - 1) * log(anew(i)) - beta * anew(i);
    double den = (alpha - 1) * log(a(i)) - beta * a(i);
    
    for(j = 0; j < Z.n_cols; j ++) {
      double plusnum = sum(b + (anew - b) % Z.col(j));
      double plusden= sum(b + (a - b) % Z.col(j));
      
      num += (lgamma(plusnum) - lgamma(sum(X.col(j)) + plusnum) + lgamma(X(i, j) + anew(i)) - lgamma(anew(i))) * Z(i, j);
      den += (lgamma(plusden) - lgamma(sum(X.col(j)) + plusden) + lgamma(X(i, j) + a(i)) - lgamma(a(i))) * Z(i, j); 
    }
    
    if(log(R::runif(0, 1)) <= (num - den)) a(i) = anew(i); else anew(i) = a(i);
  }
  
  return(a);
}

// [[Rcpp::export]]
vec updateb(mat X, mat Z, vec a, vec b, double alpha = 1, double beta = 0.1, double s = 1) {
  int i, j; vec bnew = b;
   
  for(i = 0; i < b.n_elem; i ++) {
    do {
      bnew(i) = R::rnorm(b(i), s);
    } while ((bnew(i) <= 0) | (bnew(i) >= a(i)));
    
    double num = (alpha - 1) * log(bnew(i)) - beta * bnew(i);
    double den = (alpha - 1) * log(b(i)) - beta * b(i);
    
    for(j = 0; j < Z.n_cols; j ++) {
      double plusnum = sum(bnew + (a - bnew) % Z.col(j));
      double plusden= sum(b + (a - b) % Z.col(j));
      
      num += (lgamma(plusnum) - lgamma(sum(X.col(j)) + plusnum) + lgamma(X(i, j) + bnew(i)) - lgamma(bnew(i))) * (1 - Z(i, j));
      den += (lgamma(plusden) - lgamma(sum(X.col(j)) + plusden) + lgamma(X(i, j) + b(i)) - lgamma(b(i))) * (1 - Z(i, j)); 
    }
    
    if(log(R::runif(0, 1)) <= (num - den)) b(i) = bnew(i); else bnew(i) = b(i);
  }
  
  return(b);
}