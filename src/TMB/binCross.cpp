#define TMB_LIB_INIT R_init_binCross
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  /* data section */
  DATA_MATRIX(X);
  DATA_VECTOR(y);   // binary response
  DATA_IVECTOR(g1);   //index of 1st random effect
  DATA_IVECTOR(g2);   //index of 2nd random effect
  
  
  /* Parameter section */
  PARAMETER_VECTOR(beta);
  PARAMETER(lsigmaA);
  PARAMETER(lsigmaB);
  PARAMETER_VECTOR(a); //random effects
  PARAMETER_VECTOR(b); //random effects
  
  using namespace density;

  Type aSD = exp(lsigmaA);
  Type bSD = exp(lsigmaB);

  Type nll = 0.0;     // Negative log likelihood function
  nll -= dnorm(a, Type(0), Type(1), true).sum();	// a's ~ N(0,1)
  nll -= dnorm(b, Type(0), Type(1), true).sum();	// b's ~ N(0,1)

  int n =  X.rows();
  
  vector<Type> eta = X * beta;

 
   
 for(int i = 0; i < n; i++){
   Type etai = eta(i) + a(g1(i)) * aSD +  b(g2(i)) * bSD;
   Type prob =  pnorm(etai, Type(0), Type(1));
   prob = squeeze(prob); 
   nll -= dbinom(y(i), Type(1), prob, true);
    };

  return nll;
}
