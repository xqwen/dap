#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string.h>

class parser {


 public:
  
  vector<vector<vector<double> > > geno_vec;
  vector<vector<vector<double> > > covar_vec;

  vector<vector<double> > pheno_vec;
  
  map<int, string> pheno_map;
  map<string, int> pheno_index;;
  map<int, string> geno_map;
  map<string, int> geno_rmap;

  string pheno_name;

  void process_data(char *filename);
  void output();
  
 private:
  
  void process_line(string line);

};
