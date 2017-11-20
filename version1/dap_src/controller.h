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
  
  void update();

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
  int print_all_model;

  // log10 normalizing constant
  double log10_pnorm;
  

  map<string,double> single_log10_abfv;

  // threshold 
  double snp_select_thresh;   // conditional inclusion prob.
  double size_select_thresh;  // the decay of the posteriors of a given size


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
  
  void set_print_all_model(){
    print_all_model = 1;
  }

  void set_thread(int thread){
    nthread = thread;
  }

  


 private:
  
  // initialization
  
  void load_grid(char *grid_file);
  void init();

  // computing engine
  double compute_log10_prior(map<int,int>& mcfg_map);
  size_model compute_post_model(vector<int>& best_model);
  size_model compute_post_model(int size, vector<int>& best_model);
  double conditional_est(vector<int> &control_set);
 
  void backward_check(vector<int> & best_model);


  // reporting
  void summarize_approx_posterior();


  // utility
  void parse_nmodel(Nmodel nmod);
  vector<int> get_config(int c);
  double *get_weights(vector<double>& vec);

  
};  

double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

