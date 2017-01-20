#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <map>


class SSLR_SS {
  
 private:
  

  int p; // genotype size
  int n; // sample size
  
  gsl_matrix *GtG;
  gsl_matrix *Gty;

  double yty;
    
  
  vector<double> phi2_vec; //effect size grid


 private:
  // options
  double sigma_option;  // 0 to 1, mixture fraction of mle of Sigma under the alternative model, default 0
    
  void print_matrix(gsl_matrix *M, int a, int b);
  
  
 public:
  
  // interface
  // empty constructor, assign default options
  SSLR_SS(){
    sigma_option = 0.0; 
    GtG= Gty = 0;
  }
  
  ~SSLR_SS();

  // init
  void init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_);
  
  // options
  void set_effect_vec(vector<double> &phi2);
  void set_sigma_option(double option){
    sigma_option = option;
  }

  
  double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

  double compute_log10_ABF(vector<int> &indicator);
  double compute_log10_ABF(vector<vector<int> > &indicator);
   
  

};


