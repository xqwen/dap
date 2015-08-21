#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>



class SSLR {
  
 private:
  

  int p; // genotype size
  int s; // subgroup size
  
  vector<int> nv; // sample size
  vector<int> qv; // control size
  

  

  int ep; //effective genotype size

  ///////// all of the followings just need to compute onece and use throughout
  
  vector<gsl_matrix*> Yv; // phenotype matrix
  vector<gsl_matrix*> Xgv; // genotype matrix
  vector<gsl_matrix*> Xcv; // control matrix 
  vector<gsl_matrix *> Tv;

  vector<double> hv; // inverse gamma prior 
  vector<double> mv;
  
  vector<double> sigma0_vec;

  

  ////////////  used per configuration
  gsl_matrix *eVb; 
  gsl_matrix *eVg_inv;
  gsl_matrix *Gamma;
  gsl_matrix *Sigma_inv; // Sigma^{-1} diagonal sxs matrix  
  
  
  // used per (phi, omg) value
  gsl_matrix *Wg; // effect prior  
  
  vector<double> sigma_vec;
  vector<double> omg2_vec; //effect size grid
  vector<double> phi2_vec; //effect size grid


 private:
  // options
  double sigma_option;  // 0 to 1, mixture fraction of mle of Sigma under the alternative model, default 0


 private:
  
  void compute_common();
  
  // utilites for computing residual error cov
  void compute_sigma(vector<vector<int> >& config);
  void compute_sigma_null();
  void compute_sigma_mle(vector<vector<int> >& config);
  void invert_Sigma();
  gsl_matrix *compute_residual(gsl_matrix *y, gsl_matrix *X, int size, int n,double &factor);
  
  
  // utilites for configuration specific computation
  
  void construct_Gamma(vector<vector<int> >& config, vector<int> &noz_vec);
  void construct_meta_Gamma(vector<vector<int> >& config, vector<int> &noz_vec);
 
  
  void set_Wg(double phi2, double omg2);
  
  // compute stats common for a configuration
  void compute_stats(vector<int> &noz_vec);
  
  // evaluating ABF
  double compute_log10_ABF(gsl_matrix *Wg);
  
  gsl_matrix *vec(gsl_matrix *M, int a, int b);
  gsl_matrix *kron (gsl_matrix *M, gsl_matrix *L, int a, int b);
  gsl_matrix *kron2 (gsl_matrix *M, int mr, int mc, gsl_matrix *L, int lr, int lc);

  void print_matrix(gsl_matrix *M, int a, int b);
  
  
 public:
  
  // interface
  // empty constructor, assign default options
  SSLR(){
    sigma_option = 0.0; 
  }


  // init
  void init(vector<vector<double> > & Y_in, vector<vector<vector<double> > > & Xg_in, vector<vector<vector<double> > > & Xc_in);
  
  // options
  void set_IG_prior(vector<double> & hv_in, vector<double> & mv_in);
  void set_effect_vec(vector<double> &phi2,vector<double>& omg2_vec);
  void set_sigma_option(double option){
    sigma_option = option;
  }


  vector<double> get_SNP_effect(int p);
  
  double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

  double compute_log10_ABF(vector<vector<int> > &indicator);
  
  double compute_log10_ABF(vector<vector<int> >& indicator, double phi2, double omg2);
 
  ~SSLR();
  


};


