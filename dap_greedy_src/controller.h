#include "parser.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "SSLR.h"
#include <set>
#include <map>
#include <stdio.h>

class NSNP {

 public:

  string name;
  double incl_prob;
  int cluster; //signal cluster membership
};

bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs);



class Nmodel {

 public: 
  string id;
  double prob;
  double post_score;
  int size;
};
bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs);



// set of models with same size
class size_model {

 public:
  int size;
  map<string, double> post_map;
  double log10_sum_post;
  vector<vector<int> > mvec; // models to be expanded for the next round
  void update();
  vector<int> snp_cluster;
};



class controller {
  
 private:


  FILE *outfd;
  
  parser pars;
  string gene;

  vector<double> phi2_vec;
  vector<double> omg2_vec;
  
  SSLR sslr;

  int p;
  int s;
  int q;

  vector<vector<int> > null_config;
  map<int,int> null_cfg_map;
  
  vector<vector<double> > pi_vec; //prior for configs
  
  double prior_ratio; // for estimate residual from the approximaiton
  

  int max_size; //user option, maximum model size

  
  // collection of all sizes of size_model
  vector<size_model> szm_vec; // size_model collections
  

  //candidate set/map of SNPs
  vector<int> cand_set;  // candidate SNP set for higher order models
  map<int, int> cand_map;


  // for reporting
  map<string,int> nsnp_map;
  
  vector<NSNP> nsnp_vec;
  vector<Nmodel> nmodel_vec;


  // output pip for all SNPs default only output those > 0.001
  int output_all;
  
  // log10 normalizing constant
  double log10_pnorm;
  

  map<string,double> single_log10_abfv;

  // threshold 
  double snp_select_thresh;   // conditional inclusion prob.
  double size_select_thresh;  // the decay of the posteriors of a given size

  // threshold for greedy algorithm
  int    size_limit;
  double cluster_pip_thresh; 
  int    priority_msize;
  double log10_bf_thresh;

  // ABF options
  double sslr_sigma_option;
  
  // openmp thread
  int nthread;

 public:

  int non_stop; // enumerating all models if non_stop = 1


 public:
  
  // public interface
  
  controller(char *data_file, char *grid_file);
  void run();
  void quick_est();
  void scan();
  // system setting
  
  // options
  

  void set_sslr_option(double sigma_option){
    sslr_sigma_option= sigma_option;
  }

  void set_size_limit(int size_limit_thresh){
    if(size_limit<0)
      size_limit = p;
    size_limit = size_limit_thresh;
  }

  
  void set_prior(char *prior_file);
  void set_prior(double pes=1, double lambda=0.5);

  void set_max_size(int msize){
    max_size = msize;
  }

  void set_snp_select_thresh (double thresh){
    snp_select_thresh = thresh;
  }

  void set_size_select_thresh (double thresh){
    size_select_thresh = thresh;
  }


  void set_outfile(char *filename);
    
  void set_gene(string gname){
    gene=gname;
  }

  void set_output_all(){
    output_all = 1;
  }
  
  void set_thread(int thread){
    nthread = thread;
  }

  


 private:
  
  // initialization
  
  void load_grid(char *grid_file);
  void set_default_grid();
  
  void init();

  // computing engine
  double compute_log10_prior(map<int,int>& mcfg_map);
  size_model compute_post_model_single(vector<int> &bm);
  size_model compute_post_model(int size, int use_abs_cutoff);
  double conditional_est(vector<int> &control_set);
 
  void backward_check(vector<int> & best_model);


  // reporting
  void summarize_approx_posterior();


  // utility
  void parse_nmodel(Nmodel nmod);
  vector<int> get_config(int c);
  double *get_weights(vector<double>& vec);
  double compute_average_r2 (const vector<int> & vec);
  double compute_average_r2(const vector<int> & vec1, const vector<int> & vec2);
  double compute_r2(int i, int j); 
};  

double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

