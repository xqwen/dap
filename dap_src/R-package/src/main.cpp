#include <Rcpp.h>
#include "controller.h"
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
int dap(List arg) {
  char grid_file[128];
  char data_file[128];
  char zval_file[128];
  char est_file[128];
  char ld_file[128];
  char out_file[128];
  char log_file[128];
  char gene_name[64];

  char prior_file[128];


  memset(gene_name,0,64);
  memset(grid_file,0,128);
  memset(out_file,0,128);
  memset(log_file,0,128);
  memset(data_file,0,128);
  memset(ld_file,0,128);
  memset(zval_file,0,128);
  memset(est_file,0,128);

  memset(prior_file,0,128);


  int ld_format = 1; // for correlation matrix


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

  strcpy(out_file, "output.dap");
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
    if(mystrings[i]=="d_zval")
    {
      strcpy(zval_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="ld")
    {
      strcpy(ld_file, arg[i]);
      continue;
    }
    if(mystrings[i]=="data_ld2")
    {
      strcpy(ld_file, arg[i]);
      ld_format = 2;
      continue;
    }
    if(mystrings[i]=="est")
    {
      strcpy(est_file, arg[i]);
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

    stop("Unknow option ", mystrings[i]);
  }

  controller con;


  if(strlen(data_file)!=0){
    con.initialize(data_file,grid_file);
  }else if(strlen(zval_file)!=0 && strlen(ld_file)!=0 ){
    con.initialize(zval_file, ld_file, grid_file,sample_size, ld_format);
  }else if(strlen(ld_file)!=0 && strlen(est_file)!=0 && sample_size >0 && syy>0){
    con.initialize(est_file, ld_file, grid_file, sample_size, syy,ld_format);
  }else{
    stop("Error: no suitable input data specified \n");
  }

  con.set_for_r();

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
  con.print_dap_config();

  con.run();
  return 0;

}




