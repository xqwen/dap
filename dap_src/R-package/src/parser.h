#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <gsl/gsl_matrix.h>

using namespace std;

class parser {


 public:
  
  vector<vector<vector<double> > > geno_vec;
  vector<vector<vector<double> > > covar_vec;

  vector<vector<double> > pheno_vec;
  
  map<int, string> pheno_map;
  map<string, int> pheno_index;
  map<int, string> geno_map;
  map<string, int> geno_rmap;

  string pheno_name;
  
  // summary level inof
  gsl_matrix *zval_matrix;
  gsl_matrix *ld_matrix;
  vector<double> gty_vec;
  vector<double> sxx_vec;

  void process_data(char *filename);
  void process_summary_data(char *zval_file, char *ld_file, int sample_size = -1, int ld_format=1);
  void process_summary_data2(char *effect_file, char *ld_file, int sample_size, double syy,int ld_format=1);
  void output();

  // for R use only
  void process_data(vector<double>& pheno, vector<vector<double> >& geno, string grp, string name, vector<string>& geno_name, bool regress);
  void process_summary_data2(vector<string>& snp_names, vector<double>& beta_vector, vector<double>& se_vector, vector<vector<double>>& ld, int sample_size, double syy, string name);
  void process_summary_data(vector<string>& snp_names, vector<double>& z_vector, vector<vector<double>>& ld, int sample_size, string name);



 private:
  
  void process_line(string line);
  void regress_cov(vector<double> &phenov, vector<vector<double> > &cov, vector<vector<double> > &genov );
};
