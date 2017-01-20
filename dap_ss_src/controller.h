#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "SSLR_SS.h"
#include <string>
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
  

  string gene;

  vector<double> phi2_vec;
  
  SSLR_SS sslr;

  gsl_matrix *GtG;
  gsl_matrix *Gty;
  double yty;


  int p;
  int s;
  int q;
  int n;
  
  int all_include;
  

  vector<double> beta_vec;
  vector<double> var_vec;
  vector<double> freq_vec;

  map<int, string> geno_map;
  map<string, int> geno_rmap;
  map<string, int> include_list;

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


  // input option
  int specify_yty;


  // output pip for all SNPs default only output those > 0.001
  int output_all;
  
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
  
  controller();


  void init(char *grid_file, char *data_file, char *ld_file);
  void init_prior();
  void run();
  

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
  
  void set_thread(int thread){
    nthread = thread;
  }

  void set_sample_size(int sample_size){
    n = sample_size;
  }
  
  void set_yty(double val){
    yty = val;
    specify_yty = 1;
  }

 private:
  
  // initialization
  
  void load_grid(char *grid_file);
  void load_summary(char *data_file);
  void set_inclusion_list(char* inc_file);
  void process_data();


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

