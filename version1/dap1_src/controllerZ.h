#include <string>
#include <vector>
#include <map>
#include <string>



class controllerZ {
  
 private:
  int p;
  map<string,int> geno_rmap;
  map<int, string> geno_map;

  map<string, int> use_snp; // subset filter
  vector<double> pvec; // priors	       
  vector<double> log10_BF_vec; //single snp log10_BF

  int set_filter;

 public:
  

  controllerZ(){
    set_filter = 0;
  }

  void set_subset(char *snp_list);

  void init(char *zval_file);
  void init_BF(char *bf_file);
  void set_prior(double pi1);
  void set_prior(char *prior_file);
  void run();
  



};


double log10_weighted_sum(vector<double> &vec, vector<double> &wts);
double compute_log10_BF(double z_score);
