#include <Rcpp.h>
#include "controller.h"
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

using namespace Rcpp;
using namespace std;

List summary_option_0(controller& con);
List summary_option_1(controller& con);
void print_dap_config(controller& con);

// [[Rcpp::export]]
List dap_main(int input_format, List arg, int quiet) {
  // input_format
  // 1: data_file (string)
  // 2: data x,y  (NumericMatrix, NumericVector)
  // 3: summary statistics (Vectors)
  // 4: zval (Vectors)
  
  char grid_file[128];
  char data_file[128];
  char out_file[128];
  char log_file[128];
  char gene_name[64];

  char prior_file[128];


  memset(gene_name,0,64);
  memset(grid_file,0,128);
  memset(out_file,0,128);
  memset(log_file,0,128);
  memset(data_file,0,128);

  memset(prior_file,0,128);


  // int ld_format = 1; // for correlation matrix
  // for now the R-package allows ld in format 1 only!!!

  double abf_option = -1;

  int msize =-1;

  int output_all = 0;

  double pes = 1.0;
  double pi1 = -1;
  //    double lambda = 0.5;
  double ld_control = -1;

  int size_limit = -1;

  int sample_size = -1;
  double syy = -1;

  double snp_select_thresh = -1;
  double size_select_thresh = -1;

  // alternative non-fm running options
  int run_scan = 0;
  int extract_ss = 0;

  int thread = 1;


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
    if(mystrings[i]=="n")
    {
      sample_size = arg[i];
      continue;
    }
    if(mystrings[i]=="syy")
    {
      syy = arg[i];
      continue;
    }


    // prior file
    if(mystrings[i]=="prior")
    {
      strcpy(prior_file, arg[i]);
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

    // thresholds

    if(mystrings[i]=="converg_thresh"){
      size_select_thresh = arg[i];
      continue;
    }

    if(mystrings[i]=="size_limit"){
      size_limit = arg[i];
      continue;
    }

    if(mystrings[i]=="no_size_limit"){
      size_limit = -1;
      continue;
    }


    if(mystrings[i]=="ld_control"){
      ld_control = arg[i];
      continue;
    }

    // msize option  for DAP-K
    if(mystrings[i]=="msize" || mystrings[i]=="dapk"){
      msize = arg[i];
      continue;
    }


    // openmp threads

    if(mystrings[i]=="t"){
      thread = arg[i];
      continue;
    }


    // output option
    if(mystrings[i]=="all"){
      output_all = 1;
      continue;
    }

    // gene/locus name
    if(mystrings[i]=="name"){
      strcpy(gene_name, arg[i]);
      continue;
    }

    // no finemappin, just sincle SNP analysis
    if(mystrings[i]=="scan"){
      run_scan = 1;
      continue;
    }

    if(mystrings[i]=="dump_summary"){
      extract_ss = 1;
      continue;
    }

    if(mystrings[i]=="dump_summary2"){
      extract_ss = 2;
      continue;
    }
  }

  controller con;
  
  if(input_format==1){
    con.initialize(data_file,grid_file);
  }else if(input_format==2){
    
    NumericMatrix x = arg["x"];
    NumericVector y = arg["y"];
    string name = arg["pheno_name"];
    
    int n_gene = x.ncol(), n_ind = x.nrow();
    vector<string> genonames = as<vector<string>>(colnames(x));
    vector<double> pheno = as<vector<double>>(y);
    vector<vector<double>> geno;
    geno.resize(n_gene);
    for(int i=0; i<n_gene; i++){
      NumericVector v = x(_,i);
      geno[i] = as<vector<double>>(v);
      geno[i].resize(n_ind);
    }
    con.initialize(pheno, name, geno, genonames, grid_file, true);
  }else if(input_format==3){
    
    StringVector snp_name = arg["snp_names"];
    NumericVector est     = arg["beta"];
    NumericVector est_se  = arg["se"];
    NumericMatrix ld_matrix = arg["ld"];
    string name = arg["pheno_name"];
    
    int p = snp_name.size();
    vector<string> genonames = as<vector<string>>(snp_name);
    vector<double> beta = as<vector<double>>(est);
    vector<double> se = as<vector<double>>(est_se);
    vector<vector<double>> ld;
    ld.resize(p);
    for(int i=0; i<p; i++){
      NumericVector v = ld_matrix(_,i);
      ld[i] = as<vector<double>>(v);
      ld[i].resize(p);
    }
    
    con.initialize(genonames, beta, se, ld, sample_size, syy, name);
  }else if(input_format==4){
    
    StringVector snp_name = arg["snp_names"];
    NumericVector     est = arg["zvals"];
    NumericMatrix ld_matrix = arg["ld"];
    string name = arg["pheno_name"];
    
    int p = snp_name.size();
    vector<string> genonames = as<vector<string>>(snp_name);
    vector<double> zvals = as<vector<double>>(est);
    vector<vector<double>> ld;
    ld.resize(p);
    for(int i=0; i<p; i++){
      NumericVector v = ld_matrix(_,i);
      ld[i] = as<vector<double>>(v);
      ld[i].resize(p);
    }
    
    con.initialize(genonames, zvals, ld, sample_size, name);
  }

  con.set_for_r(quiet);

  con.set_outfile(out_file, log_file);

  con.set_gene(gene_name);
  con.set_abf_option(abf_option);
  con.set_thread(thread);

  con.set_size_limit(size_limit);

  if(ld_control>=0)
    con.set_ld_control(ld_control);

  if(msize>=1){
    con.set_max_size(msize);
  }


  if(strlen(prior_file)==0){
    if(pi1 != -1){
      if(0<pi1 && pi1 <1){
        con.set_prior(pi1);
      }else{
        warning("Warning: pi1 specification is outside the range, ignored...\n");
        con.set_prior_exp(pes);
      }
    }else{
      // default
      con.set_prior_exp(pes);
    }
  }else
    con.set_prior(prior_file);

  if(output_all == 1)
    con.set_output_all();

  if(snp_select_thresh>=0)
    con.set_snp_select_thresh(snp_select_thresh);

  if(size_select_thresh >=0)
    con.set_size_select_thresh(size_select_thresh);

  con.run_option = 0;

  if(run_scan){
    con.run_option = 1;
  }

  if(extract_ss==1){
    con.run_option =2;
  }
  if(extract_ss==2){
    con.run_option =3;
  }


  // all done, print all configs
  if(quiet==0) print_dap_config(con);

  con.run();

  if(con.run_option==0) return summary_option_0(con);
  if(con.run_option==1) return summary_option_1(con);

  // API not built for other options
  return List::create();

}


List summary_option_0(controller& con){
  int n_nmodel = con.nmodel_vec.size();
  NumericVector   nmodel_prob, nmodel_post;
  IntegerVector   nmodel_size;
  CharacterVector nmodel_name;

  for(int i=0;i<n_nmodel;i++){
    Nmodel& nm = con.nmodel_vec[i];

    nmodel_prob.push_back(nm.prob);
    nmodel_size.push_back(nm.size);
    nmodel_post.push_back(nm.post_score);
    nmodel_name.push_back(nm.id);
  }

  DataFrame nmodel = DataFrame::create(Named("posterior")     = nmodel_prob,
                                       Named("size")          = nmodel_size,
                                       Named("score")         = nmodel_post,
                                       Named("configuration") = nmodel_name);
  NumericVector msize = NumericVector::create(Named("mean")=con.get_msize_mean(), Named("sd")=con.get_msize_sd());


  double min_pip    = con.get_min_pip();
  int    output_all = con.get_output_all();
  int    n_snp      = con.nsnp_vec.size();

  CharacterVector nsnp_name;
  NumericVector   nsnp_prob, nsnp_abfv;
  IntegerVector   nsnp_cluster;

  for(int i=0;i<n_snp;i++){
    NSNP& nsnp = con.nsnp_vec[i];

    if(nsnp.incl_prob < min_pip)
      nsnp.incl_prob = min_pip;
    if(nsnp.incl_prob<1e-3&&!output_all)
      break;
    if(nsnp.cluster==-1&&!output_all)
      continue;
    nsnp_name.push_back(nsnp.name);
    nsnp_prob.push_back(nsnp.incl_prob);
    nsnp_abfv.push_back(con.single_log10_abfv[nsnp.name]);
    nsnp_cluster.push_back(nsnp.cluster);
  }

  DataFrame SNP = DataFrame::create(Named("predictor")=nsnp_name,
                                    Named("pip")=nsnp_prob,
                                    Named("log10abf")=nsnp_abfv,
                                    Named("cluster") =nsnp_cluster);

  List model_summary = List::create(Named("model") = nmodel,
                                    Named("model.size") = msize,
                                    Named("log10NC") = con.get_log10_pnorm(),
                                    Named("PIP.min") = min_pip,
                                    Named("N") = con.get_N(),
                                    Named("response") = con.get_pheno_name());


  List result = List::create(Named("variant") = SNP,
                             Named("model.summary") = model_summary);

  if(con.grp_vec.size() > 0){
    DataFrame cluster = DataFrame::create(Named("cluster.size") = wrap(con.cluster_count),
                                          Named("cluster.pip") = wrap(con.cluster_pip),
                                          Named("average.r2")  = wrap(con.cluster_r2));

    int n_cluster = con.cluster_count.size();
    NumericMatrix r2 (n_cluster);
    for(int i=0; i<n_cluster; i++){
      r2(i,i) = con.cluster_r2[i];
      int counti = con.grpr2_map[i];
      for(int k=i+1; k<n_cluster; k++){
        int countk = con.grpr2_map[k];
        string id;
        if(countk<counti)
          id = to_string(countk)+ ":"+to_string(counti);
        if(countk>counti)
          id = to_string(counti)+":"+to_string(countk);
        r2(i,k) = con.grp_r2[id];
        r2(k,i) = r2(i,k);
      }
    }
    rownames(r2) = wrap(con.cluster_id);
    colnames(r2) = rownames(r2);

    List signal_cluster = List::create(Named("cluster.summary") = cluster,
                                       Named("cluster.r2") = r2);

    result.push_front(signal_cluster, "signal.cluster");
  }

  return result;
}

List summary_option_1(controller& con){
  DataFrame result = DataFrame::create(Named("SNP")=con.single_snp_name,
                                       Named("log10_ABF")=con.single_snp_log10_ABF,
                                       Named("b")=con.single_snp_bhat,
                                       Named("se")=con.single_snp_se,
                                       Named("z-value")=con.single_snp_zval);
  return List::create(result);
}

void print_dap_config(controller& con){
  Rcout << "\n============ DAP Configuration ============\n\n";
  Rcout << "INPUT\n\n";

  Rcout << "\t* individual-level data\n";
  Rprintf("\t* number of candidate predictors: %i\n", con.get_p());
  Rprintf("\t* sample size: %i\n\n", con.get_N());

  Rcout << "PROGRAM OPTIONS\n\n";
  Rprintf("\t* maximum model size allowed [ msize = %i ] ", con.get_max_size());
  if(con.get_max_size()==con.get_p())
    Rcout << "(no restriction)\n";
  else
    Rcout << "\n";

  Rprintf("\t* LD control threshold [ ld_control = %.2f ]\n", con.get_ld_control_thresh());
  Rprintf("\t* normalizing constant convergence threshold [ converg_thresh = %.2e ] (log10 scale)\n", con.get_size_select_thresh());

  #ifdef _OPENMP
    Rprintf("\t* number of parallel threads [ thread = %i ]\n", con.get_nthread());
  #else
    Rcout << "\t* number of parallel threads [ thread = 1 ] (OpenMP not available)\n";
  #endif

  Rcout << "\n===========================================\n";
}


