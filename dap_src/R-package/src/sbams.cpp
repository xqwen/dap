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


// [[Rcpp::export]]
List extract_sbams(List arg){
  char grid_file[128];
  char data_file[128];
  char out_file[128];
  char log_file[128];
  char gene_name[64];
  
  
  memset(gene_name,0,64);
  memset(grid_file,0,128);
  memset(out_file,0,128);
  memset(log_file,0,128);
  memset(data_file,0,128);
  
  double abf_option = -1;

  bool has_prior_data = false;
  double pes = 1.0;
  double pi1 = -1;
  
  
  vector< string > mystrings =  arg.attr("names");
  
  for(int i=0; i<arg.size(); i++)
  {
    // required data files and additional info
    if(mystrings[i]=="data")
    {
      strcpy(data_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="grid")
    {
      strcpy(grid_file, arg[i]);
      continue;
    }
    
    
    // prior file
    if(mystrings[i]=="has_prior"){
      has_prior_data = true;
      continue;
    }
    
    // prior options
    
    if(mystrings[i]=="ens"){
      pes = arg[i];
      continue;
    }
    
    
    if(mystrings[i]=="pi1"){
      pi1 = arg[i];
      continue;
    }
    
    if(mystrings[i]=="abf"){
      abf_option = arg[i];
      continue;
    }
    
    // gene/locus name
    if(mystrings[i]=="name"){
      strcpy(gene_name, arg[i]);
      continue;
    }
    
  }
  
  controller con;
  con.initialize(data_file,grid_file);
  
  con.set_for_r(1);
  
  con.set_outfile(out_file, log_file);
  
  con.set_gene(gene_name);
  con.set_abf_option(abf_option);


  if(pes>0){
    con.set_ens(pes);
  }
  
  
  if(!has_prior_data){
    if(pi1 != -1){
      if(0<pi1 && pi1 <1){
        con.set_prior(pi1);
      }else{
        // fprintf(stderr, "Warning: pi1 specification is outside the range, ignored...\n");
        con.set_prior_default();
      }
    }else{
      // default
      con.set_prior_default();
    }
  }else{
    // StringVector  prior_snp_names_r = arg["prior_snp_names"];
    // NumericVector prior_values_r    = arg["prior_values"];
    vector<string> prior_snp_names = as<vector<string>>(arg["prior_snp_names"]);
    vector<double> prior_values    = as<vector<double>>(arg["prior_values"]);
    con.set_prior(prior_snp_names, prior_values);
  }
  
  con.extract_ss2_in_r();
  
  int p = con.get_p();
  NumericVector ld_matrix_flat = wrap(con.ld_matrix_flat);
  ld_matrix_flat.attr("dim") = Dimension(p,p);
  NumericMatrix ld_matrix = as<NumericMatrix>(ld_matrix_flat);
  
  DataFrame est = DataFrame::create(Named("snp")=wrap(con.single_snp_name),
                                    Named("b")  =wrap(con.single_snp_bhat),
                                    Named("se") =wrap(con.single_snp_se));
  
  return List::create(Named("est")=est,
                      Named("ld") =ld_matrix,
                      Named("n")  =con.get_N(),
                      Named("syy")=con.get_syy(),
                      Named("name")=con.get_pheno_name());
}



