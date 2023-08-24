#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>


// [[Rcpp::export]]
double likAGH(double rho, Rcpp::List list_eta, Rcpp::List list_w, int niter, DoubleVector ws, DoubleVector z)
{
  int len = list_w.size();
  double sigma2 = rho / (1.0 - rho);
  double sigma = sqrt(sigma2);
  double ll = 0.0;
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    IntegerVector wi = list_w[i];
    int ni = wi.size();
    double ui = 0.0;
    double H;
    for (int k=0; k<niter; k++){
      double g = -1.0 * ui;
      H = 1.0;
      for (int j=0; j<ni; j++){
        double arg = wi[j] * (etai[j] * sqrt(1.0 + sigma2) + sigma * ui);
        double phi = R::dnorm(arg, 0.0, 1.0, 0);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
        g += wi[j] * sigma * phi / Phi;
        H += pow(wi[j], 2.0) * sigma2 * phi / pow(Phi, 2.0) * (arg  * Phi + phi);
      }
      ui += g / H;
    }
   double se = 1.0 / sqrt(H);
   int nq = ws.size();
   // now computes AGH
   double Ii = 0.0;
   for (int k=0; k<nq; k++){
   double zk = sqrt(2) * se * z[k] + ui;
   double f = 0.0;
   for (int j=0; j<ni; j++){
      double arg = wi[j] * (etai[j] * sqrt(1.0 + sigma2) + sigma * zk);
      f += R::pnorm(arg, 0.0, 1.0, 1, 1);
      }
    f += R::dnorm(zk, 0.0, 1.0, 1);
    Ii += ws[k] * exp(f);
     }
    Ii *= se *   sqrt(2.0);
    //If Ii == 0 then switch to Laplace
    if(Ii == 0) {
      double zk = ui;
      double f = 0.0;
      for (int j=0; j<ni; j++){
        double arg = wi[j] * (etai[j] * sqrt(1.0 + sigma2) + sigma * zk);
        f += R::pnorm(arg, 0.0, 1.0, 1, 1);
       }
      f += R::dnorm(zk, 0.0, 1.0, 1);
      ll += log(1.772454) + f + log(se) + log(sqrt(2.0));
     }
    else ll += log(Ii);
   }
  return(ll);
}



// [[Rcpp::export]]
DoubleVector getEffects(double rho, Rcpp::List list_eta, Rcpp::List list_w, int niter, DoubleVector ws, DoubleVector z)
{
  int len = list_w.size();
  double sigma2 = rho / (1.0 - rho);
  double sigma = sqrt(sigma2);
  DoubleVector  out (len);
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    IntegerVector wi = list_w[i];
    int ni = wi.size();
    double H = 0.0;
    double ui = 0.0;
    for (int k=0; k<niter; k++){
      double g = -1.0 * ui;
      H = 1.0;
      for (int j=0; j<ni; j++){
        double arg = wi[j] * (etai[j] * sqrt(1.0 + sigma2) + sigma * ui);
        double phi = R::dnorm(arg, 0.0, 1.0, 0);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
        g += wi[j] * sigma * phi / Phi;
        H += pow(wi[j], 2.0) * sigma2 * phi / pow(Phi, 2.0) * (arg  * Phi + phi);
      }
      ui += g / H;
    }
     out[i] = ui;
  }
  return(out);
}


// [[Rcpp::export]]
DoubleVector getEffectsOrd(double rho, DoubleVector alphae, Rcpp::List list_eta, Rcpp::List list_y, int niter, DoubleVector ws, DoubleVector z)
{
  int len = list_y.size();
  double sigma2 = rho / (1.0 - rho);
  double sigma = sqrt(sigma2);
  DoubleVector  out (len);
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    IntegerVector yi = list_y[i];
    int ni = yi.size();
    double H = 0.0;
    double ui = 0.0;
    for (int k=0; k<niter; k++){
      double g = -1.0 * ui;
      H = 1.0;
      for (int j=0; j<ni; j++){
	double arg = (alphae[yi[j]] - etai[j]) * sqrt(1.0 + sigma2) - sigma * ui;
	double arg1 = (alphae[yi[j] - 1] - etai[j]) * sqrt(1.0 + sigma2) - sigma * ui;
        double phi = R::dnorm(arg, 0.0, 1.0, 0);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
	double phi1 = R::dnorm(arg1, 0.0, 1.0, 0);
        double Phi1 = R::pnorm(arg1, 0.0, 1.0, 1, 0);
	g -= sigma * (phi - phi1) / (Phi - Phi1);
        H += sigma2 * (((arg * phi - arg1 * phi1 ) * (Phi - Phi1)) + pow(phi - phi1, 2.0)) / pow(Phi - Phi1, 2.0);
      }
     ui += g / H;
    }
     out[i] = ui;
  }
  return(out);
}



// [[Rcpp::export]]
double likAGHOrd(double rho,  DoubleVector alphae, Rcpp::List list_eta, Rcpp::List list_y, int niter, DoubleVector ws, DoubleVector z)
{
  int len = list_y.size();
  double sigma2 = rho / (1.0 - rho);
  double sigma = sqrt(sigma2);
  double ll = 0.0;
  // first locate the modes in u
  for(int i=0; i<len; i++){
    DoubleVector etai = list_eta[i];
    IntegerVector yi = list_y[i];
    int ni = yi.size();
    double H = 0.0;
    double ui = 0.0;
    for (int k=0; k<niter; k++){
      double g = -1.0 * ui;
      H = 1.0;
      for (int j=0; j<ni; j++){
	double arg = (alphae[yi[j]] - etai[j]) * sqrt(1.0 + sigma2) - sigma * ui;
        double arg1 = (alphae[yi[j] - 1] - etai[j]) * sqrt(1.0 + sigma2) - sigma * ui;
        double phi = R::dnorm(arg, 0.0, 1.0, 0);
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
	double phi1 = R::dnorm(arg1, 0.0, 1.0, 0);
        double Phi1 = R::pnorm(arg1, 0.0, 1.0, 1, 0);
	g -= sigma * (phi - phi1) / (Phi - Phi1);
        H += sigma2 * (((arg * phi - arg1 * phi1 ) * (Phi - Phi1)) + pow(phi - phi1, 2.0)) / pow(Phi - Phi1, 2.0);
      }
     ui += g / H;
    }
   double se = 1.0 / sqrt(H);
   int nq = ws.size();
   // now computes AGH
   double Ii = 0.0;
   for (int k=0; k<nq; k++){
   double zk = sqrt(2) * se * z[k] + ui;
   double f = 0.0;
   for (int j=0; j<ni; j++){
     double arg  = (alphae[yi[j]] - etai[j]) * sqrt(1.0 + sigma2) - sigma * zk;
     double arg1 = (alphae[yi[j] - 1] - etai[j]) * sqrt(1.0 + sigma2) - sigma * zk;
      double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
      double Phi1 = R::pnorm(arg1, 0.0, 1.0, 1, 0);
      f += log(Phi - Phi1);
      }
    f += R::dnorm(zk, 0.0, 1.0, 1);
    Ii += ws[k] * exp(f);
     }
    Ii *= se * sqrt(2.0);
    //If Ii == 0 then switch to Laplace
    if(Ii == 0) {
      double zk = ui;
      double f = 0.0;
      for (int j=0; j<ni; j++){
        double arg  = (alphae[yi[j]] - etai[j]) * sqrt(1.0 + sigma2) - sigma * zk;
        double arg1 = (alphae[yi[j] - 1] - etai[j]) * sqrt(1.0 + sigma2) - sigma * zk;
        double Phi = R::pnorm(arg, 0.0, 1.0, 1, 0);
        double Phi1 = R::pnorm(arg1, 0.0, 1.0, 1, 0);
        f += log(Phi - Phi1);
       }
      f += R::dnorm(zk, 0.0, 1.0, 1);
      ll += log(1.772454) + f + log(se) + log(sqrt(2.0));
     }
     else ll += log(Ii);
     }
  return(ll);
}






