#ifndef __CLASSDEF_H_
#define __CLASSDEF_H_

#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


using namespace std;




class SNP {

 public:

  string id;
  int index;
  double log10_BF;


  SNP(string snp_id, double snp_log10_BF, int snp_index){
    id = snp_id;
    log10_BF = snp_log10_BF;
    index = snp_index;
  }
 
    
 
};





class Locus {
  
 public:
  
  // possible gene information/annotation
  string id;
  
  vector<SNP> snpVec;
  
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;
  
  double log10_lik; // log10 of marginal likelihood
  
  double fdr;
  
  
  Locus(string locus_id,  vector<SNP> & snpVec_){ id = locus_id;  snpVec = snpVec_; prior_vec = pip_vec =0;  };
  Locus(){};

  
  void EM_update();
  void compute_fdr();


};






class controller {

 public:
 
  controller(){
    p=kc=kd=dist_bin_level = 0;
    force_logistic = 0;
    dist_bin_size = -1;
    finish_em = 0;
  }

  // storage
  vector<Locus> locVec;
  
  
  int p; // number of loc-SNP pairs
  
  int kc; // number of continuous covariate
  int kd; // number of discrete covariate


  double dist_bin_size;


  map<string,int> loc_hash;
  map<string,int> snp_hash;
  vector<string> cvar_name_vec;
  vector<string> dvar_name_vec;



  
  gsl_vector_int *dist_bin;
  map<int, int> dtss_map;
  map<int, int> dtss_rmap;
  int dist_bin_level;
  double EM_thresh;


  
  gsl_vector_int *dlevel; // (kd+1) entry levels of each factor


  gsl_matrix *Xc;  // p x kc
  gsl_matrix_int *Xd; // p x kd
  
  
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;
  gsl_vector *beta_vec;
  
  int ncoef;

  double final_log10_lik;
  


  int nthread; 


  int finish_em;


  void load_data(char *filename);    // load data with MatrixeQTL format -- default
  void load_data_BF(char *filename); // load data with pre-computed log10 Bayes factors

  void load_map(char *gene_map_file, char *snp_map_file);
  void load_annotation(char *annot_file);
  int count_factor_level(int col);
  
  void simple_regression();
  void single_ct_regression();
  int force_logistic;
  

  void init_params();
    
  void run_EM();

  
  void find_eGene(double thresh=0.05);
  void estimate();
  void dump_prior(char *path);


};




double log10_weighted_sum(vector<double> & val_vec, vector<double> & wts_vec);
double compute_log10_BF(double beta, double se_beta);
bool   rank_by_fdr (const Locus & loc1 , const Locus & loc2);
  
int classify_dist_bin(int snp_pos, int tss, double bin_size = -1);

#endif
