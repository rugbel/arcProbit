#define TMB_LIB_INIT R_init_ordCross
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* data section */
  DATA_MATRIX(X);
  DATA_IVECTOR(y);   // response, integer from 1 to d 
  DATA_INTEGER(d);   //number of categories
  DATA_SCALAR(maxtau); //constant to use as max tau
  DATA_IVECTOR(g1);   //index of 1st random effect
  DATA_IVECTOR(g2);   //index of 2nd random effect
  
  /* Parameter section */
  PARAMETER_VECTOR(tau);
  PARAMETER_VECTOR(beta);
  PARAMETER(lsigmaA);
  PARAMETER(lsigmaB);
  PARAMETER_VECTOR(a); //random effects
  PARAMETER_VECTOR(b); //random effects
 
  
  using namespace density;

  Type aSD = exp(lsigmaA);
  Type bSD = exp(lsigmaB);
  
  
  int n = X.rows();

  Type nll = 0.0;     // Negative log likelihood function
  nll -= dnorm(a, Type(0), Type(1), true).sum();	// a's ~ N(0,1)
  nll -= dnorm(b, Type(0), Type(1), true).sum();	// b's ~ N(0,1)
  

  vector<Type>  zetae(d+1);
  zetae(0) = -maxtau;
  for(int j=1;j<d;j++)
      zetae(j) = tau(j-1);
   zetae(d) = maxtau; 
   
 vector<Type> eta = X * beta;
   
   
  for(int i=0;i<n;i++){
    Type etai = eta(i) + a(g1(i)) * aSD +  b(g2(i)) * bSD; 
    Type up = zetae(y(i)) - etai;
    Type low = zetae(y(i) - 1) - etai;
    Type prob =  pnorm(up,Type(0),Type(1)) - pnorm(low,Type(0),Type(1));
    prob = squeeze(prob); 
    nll -= log(prob);
   };

  return nll;
}
