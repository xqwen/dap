#include <Rcpp.h>
#include "controller.h"
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cstring>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame read_sbams(const char * data_file){
  char filename[128];
  memset(filename,0,128);
  strcpy(filename, data_file);

  parser pars;
  pars.process_data(filename);

  if(pars.pheno_vec.size()==0)
    stop("[Error] No phenotype information found in the file!\n");
  if(pars.geno_vec.size()==0)
    stop("[Error] No genotype information found in the file!\n");

  if(pars.pheno_vec.size()!=1 || pars.geno_vec.size()!=1)
    stop("Please contact the author for support of multi-group sbams files.\n");

  List result_list = wrap(pars.geno_vec[0]);
  result_list.names() = wrap(pars.geno_map);
  result_list.push_front(wrap(pars.pheno_vec[0]), pars.pheno_name);

  DataFrame result = result_list;
  return result;
}



