// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]
#include <mvtnormAPI.h>


// [[Rcpp::export]]
double pmvnorm_cpp(arma::vec& upper, double rho, double abseps = 1e-3)
{
  int n = 2;
  int nu = 0;
  int maxpts = 25000;      // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  double* upper_ = upper.memptr();
  int infin[2] = {0,0};
  double lower[2] = {0,1};
  double delta[2] = {0,0};
  double corr[1];
  corr[0] = rho;
  double error = 0;
  double value = 0;
  int inform = 0;
  
  mvtnorm_C_mvtdst(&n, &nu, lower, upper_, infin, corr, delta, &maxpts, 
                   &abseps, &releps, &error, &value, &inform, &rnd);
  return(value);
}



// [[Rcpp::export]]
double lik2_binary(double rho, Rcpp::NumericVector etaG, Rcpp::List list_ind, Rcpp::IntegerVector y)
{
  int len = list_ind.size();
  double nll = 0;
  
  for(int i=0; i<len; i++){
    Rcpp::IntegerVector  veci = list_ind[i]; 
    int ni = veci.size();
    for(int j = 0; j < ni; j++)
    {
      int indfirst = veci(j) - 1;
      int y1 = y(indfirst); 
      for(int h = j+1; h < ni; h++)
      {
        int indsecond = veci(h) - 1;
        int y2 = y(indsecond);
        arma::vec upper = {etaG(indfirst), etaG(indsecond)};
        double I12 = pmvnorm_cpp(upper, rho);
        if(y1 == 1 & y2 == 1) nll += log(I12);
        if(y1 == 1 & y2 == 0) nll += log(R::pnorm(etaG(indfirst), 0.0, 1.0, 1, 0) - I12);
        if(y1 == 0 & y2 == 1) nll += log(R::pnorm(etaG(indsecond), 0.0, 1.0, 1, 0) - I12);
        if(y1 == 0 & y2 == 0) nll += log(1- R::pnorm(etaG(indfirst), 0.0, 1.0, 1, 0) -  R::pnorm(etaG(indsecond), 0.0, 1.0, 1, 0) + I12);
      }
    }
  }
  return(nll);
}



// [[Rcpp::export]]
double lik2_ord(double rho, Rcpp::NumericVector tauG, Rcpp::NumericVector etaG, Rcpp::List list_ind, 
                Rcpp::IntegerVector y)
{
  int len = list_ind.size();
  double nll = 0;
  
  for(int i=0; i<len; i++){
    Rcpp::IntegerVector  veci = list_ind[i]; 
    int ni = veci.size();
    for(int j = 0; j < ni; j++)
    {
      int indfirst = veci(j) - 1;
      int y1 = y(indfirst); 
      double tau1 = tauG(y1) - etaG(indfirst);
      double tau1m = tauG(y1 - 1) - etaG(indfirst);
      for(int h = j+1; h < ni; h++)
      {
        int indsecond = veci(h) - 1;
        int y2 = y(indsecond);
        double tau2 = tauG(y2) - etaG(indsecond);
        double tau2m = tauG(y2-1) - etaG(indsecond);
        arma::vec upper11 = {tau1m, tau2m};
        double I11 = pmvnorm_cpp(upper11, rho);
        arma::vec upper12 = {tau1m, tau2};
        double I12 = pmvnorm_cpp(upper12, rho);
        arma::vec upper21 = {tau1, tau2m};
        double I21 = pmvnorm_cpp(upper21, rho);
        arma::vec upper22 = {tau1, tau2};
        double I22 = pmvnorm_cpp(upper22, rho);
        nll += log(I22 + I11 - I21 - I12);
      }
    }
  }
  return(nll);
}



