using namespace std;

#include "controller.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>

#include <omp.h>

controller::controller(char *data_file, char *grid_file){

  load_grid(grid_file);

  pars.process_data(data_file);
 
  
  p = pars.geno_vec[0].size();
  s = pars.pheno_vec.size();
  q = pars.covar_vec[0].size();
   
  
  // default maximum model size
  max_size = p;
  // default thread 
  nthread = 1;

  output_all = 0;
  size_select_thresh = 0.01;  // set to big positive numbers to enforce the adaptive stopping rule
  snp_select_thresh = 0.01;


}


void controller::load_grid(char *grid_file){

  ifstream gfile(grid_file);
  string line;
  istringstream ins;

  while(getline(gfile,line)){

    ins.clear();
    ins.str(line);
    double phi;
    double omg;
    if(ins>>phi>>omg){
      phi2_vec.push_back(pow(phi,2));
      omg2_vec.push_back(pow(omg,2));
    }
  }
  gfile.close();

}



void controller::set_outfile(char *outfile){
  
  if(strlen(outfile)==0)
    outfd = stdout;
  else
    outfd = fopen(outfile,"w");
}

void controller::set_prior(double pes, double lambda){
  
  int size = 1<<s;
  
  // no configuration option yet!!!!!
  lambda = 1;
  vector<double> pvec;
  double pi1 = pes/double(p);
  if(pi1>=1)
    pi1 = 0.5;
  for(int i=0;i<size;i++){
    if(i==0){
      pvec.push_back(1-pi1); 
    }else if(i==size-1){
      pvec.push_back(pi1*lambda);
    }else{
	pvec.push_back(pi1*(1-lambda)/double(size-2));
    }
  }
  
  for(int i=0;i<p;i++){
    pi_vec.push_back(pvec);
  }
  
  return;
}

void controller::set_prior(char *prior_file){

  if(prior_file==0)
    return;

  int size = 1<<s;

  ifstream ifile(prior_file);
  string line;
  istringstream ins;

  for(int i=0;i<p;i++){
    vector<double> pvec(size,0.0);
    pi_vec.push_back(pvec);
  }

  while(getline(ifile,line)){

    ins.clear();
    ins.str(line);
    double value;
    string snp_name;
    ins>>snp_name;

    int index = pars.geno_rmap[snp_name];
    ins>>value;
    pi_vec[index][0] = 1-value;
    pi_vec[index][size-1] = value;

  }

  return;

  
}



void controller::init(){  
  
  vector<vector<vector<double> > > Xcv;
  vector<vector<vector<double> > > Xgv;

   
  Xgv = pars.geno_vec;
  Xcv = pars.covar_vec;
      
  if(Xgv.size() == 0)
    exit(0);
  
  
  // set SSLR parameters and options
  sslr.set_sigma_option(sslr_sigma_option); 
  sslr.init(pars.pheno_vec,Xgv,Xcv);
  
  //sslr.set_prior_option(sslr_prior_option);
  sslr.set_effect_vec(phi2_vec,omg2_vec);
  

  vector<int> vec0(s,0);
   
  //vector<int> cfg0;
  int total_cfg = (1<<s);
  
  //start with empty table
  for(int i=0;i<p;i++){
    null_config.push_back(vec0);
    null_cfg_map[i] = 0;
  }
  
  double sum = 0;
  for(int i=0;i<p;i++){
    sum += (1-pi_vec[i][0])/pi_vec[i][0];
  }
  
  prior_ratio = sum/p;

}



 
double controller::compute_log10_prior(map<int,int> & mcfg_map){
  
  double lp=0;
  for(int i=0;i<p;i++){
    lp += log(pi_vec[i][mcfg_map[i]]);
  }
  
  return lp/log(10);

}


void controller::quick_est(){
}


void controller::run(){
  

  init();
  
  vector<double> log10_pmass_vec;



  double curr_val = compute_log10_prior(null_cfg_map);
  log10_pmass_vec.push_back(curr_val);

  fprintf(stderr, "\n\n\tModel_Size \tcandidates \t  log10(NC)\n");

  // record best model
  vector<int> bm;
  // single SNP scan, label candidate SNPs for higher order models
  szm_vec.push_back(compute_post_model(bm));

  log10_pmass_vec.push_back(szm_vec[0].log10_sum_post);


  vector<double> wv(log10_pmass_vec.size(),1.0);
  double val = log10_weighted_sum(log10_pmass_vec,wv);
  fprintf(stderr,  "\t %4d \t \t    -\t \t %9.3f\n",1,val);


  // start exploring higher-order model
  int cs = 2;
  double prev_val = val;
  while(1){

    // user specified stopping K (DAP-K)
    if(cs>=max_size){
      fprintf(outfd, "\n*** Warning: larger models *may* contain substantial probability mass, increase maximum model size and re-run ***\n");
      break;
    }

    fprintf(stderr,"\t %4d \t \t ",cs);
    fflush(stderr);

    // finding candidate SNP sets
    conditional_est(bm);


    if(cand_set.size()<cs)
      break;
    bm.clear();
    
    
    double cps = szm_vec[szm_vec.size()-1].log10_sum_post;
    // core computation for model size = cs
    szm_vec.push_back(compute_post_model(cs,bm));


    log10_pmass_vec.push_back(szm_vec[szm_vec.size()-1].log10_sum_post);


    // update log10(NC)
    vector<double> nwv(log10_pmass_vec.size(),1.0);
    val = log10_weighted_sum(log10_pmass_vec,nwv);
    fprintf(stderr,  "%4d \t \t %9.3f\n",int(cand_set.size()),val);


    // stopping criteria
    double ncps = szm_vec[szm_vec.size()-1].log10_sum_post;
    double rb = log10(double(p)-cs+1)+log10(prior_ratio) + cps;
    double lb = log10(double(p-2*cs+2)/cs) + log10(prior_ratio) + cps;
    if( ncps <= rb && val - prev_val <= size_select_thresh){
      backward_check(bm);
      break;
    }

    prev_val = val;
    cs++;

  }


  fprintf(stderr, "\n\n");
  
  summarize_approx_posterior();
  
}



  
double controller::conditional_est(vector<int>& bm){
  
  /*
  printf("current best model: [");
  for(int i=0;i<bm.size();i++){
    printf("%d::",bm[i]);
  }
  printf("]\n");
  */
  
    
  int totalc = (1<<s)-1;
  double max = -9999;
  
  vector<vector<int> > mcfg = null_config;
  map<int,int> mcfg_map = null_cfg_map;
  
  map<int, int> control_map;
  for(int i=0;i<bm.size();i++){
    int index = bm[i];
    mcfg[index] = get_config(totalc);
    mcfg_map[index] = totalc;
    control_map[index] = 100;
  
  }
  
  double log10_baseline = sslr.compute_log10_ABF(mcfg);

  vector<double> rstv;
  vector<double> abf_vec;
  for(int i=0;i<p;i++){

    if(control_map[i]==100){
      abf_vec.push_back(0);
      rstv.push_back(log10_baseline + compute_log10_prior(mcfg_map));
      continue;
    }

    vector<vector<int> > cmcfg = mcfg;
    map<int,int> cmcfg_map = mcfg_map;
   
        
    cmcfg[i] = get_config(totalc);
    cmcfg_map[i] = totalc;
    double log10_abf =  sslr.compute_log10_ABF(cmcfg);
    
    
    abf_vec.push_back(log10_abf-log10_baseline);
    double rst = log10_abf + compute_log10_prior(cmcfg_map);
    //printf("add %d => %.3f \n", i, rst);
    rstv.push_back(rst);
    if(rst > max){
      max=rst;
    }
  }

  

  vector<double> wv(p,1.0);
  double log10_sum = sslr.log10_weighted_sum(abf_vec,wv);
  

  vector<double> post_vec;
  
  double KL = 0;
  double maxp = 0;
  int max_id = -1;
  int csize = cand_map.size();
  for(int i=0;i<p;i++){
    double val = pow(10, abf_vec[i]-log10_sum);
    if(val > maxp){
      maxp = val;
      max_id = i;
    }
    KL += (1.0/p)*log(1.0/(p*val));
    if(val>snp_select_thresh){
      cand_map[i] = 1;
      post_vec.push_back(rstv[i]);
    }
  }
  
  // if nothing added, add the best SNP
  
  if(cand_map.size()==csize&&bm.size()==1){
    cand_map[max_id] = 1;
  }
  
  //printf("KL = %f \n",KL);
  
 
  if(post_vec.size()==0)
    post_vec.push_back(max);
  
  vector<double> pwv(post_vec.size(),1.0);
  

  map<int, int>::iterator iter;
  cand_set.clear();
  for(iter=cand_map.begin();iter!=cand_map.end();iter++){
    cand_set.push_back(iter->first);
  }
  
  return sslr.log10_weighted_sum(post_vec,pwv);
}
  


// single SNP model, required to run 
size_model controller::compute_post_model(vector<int>& best_model){
  
  int totalc = (1<<s)-1;

  size_model smod;
  smod.size = 1;
  smod.log10_sum_post = 0;
  vector<double> post_vec;
  vector<double> abf_vec;
  
  double max_log10_post = -9999;

  for(int index=0;index<p;index++){
    
    string name = pars.geno_map[index];
    vector<vector<int> > mcfg = null_config;
    map<int,int> mcfg_map = null_cfg_map;
    mcfg[index]=get_config(totalc);
    mcfg_map[index] = totalc;
    double log10_abf = sslr.compute_log10_ABF(mcfg); 
    single_log10_abfv[name] = log10_abf;
    abf_vec.push_back(log10_abf);
    double log10_post =  log10_abf + compute_log10_prior(mcfg_map);
    post_vec.push_back(log10_post);
    
    
    if(log10_post>max_log10_post){
      max_log10_post = log10_post;
      best_model.clear();
      best_model.push_back(index);
    }
    smod.post_map[name] = log10_post;
    
  }

  
  
  vector<double> wv(post_vec.size(),1.0);
  smod.log10_sum_post = sslr.log10_weighted_sum(post_vec,wv);
  
  double log10_sum = sslr.log10_weighted_sum(abf_vec,wv);

  double maxp = 0;
  int max_id = -1;
  int csize = cand_set.size();
  //printf("1 cand_set: ");
  for(int i=0;i<p;i++){
    double val = pow(10, abf_vec[i]-log10_sum);
    if(val>maxp){
      maxp = val;
      max_id = i;
    }
    if(val>snp_select_thresh){
      cand_map[i] = 1;
      //printf("%d ",i);
    }
  }
  //if nothing added, add the top SNP
  if(cand_map.size()==csize){
    cand_map[max_id] = 1;
  }

  //printf("\n");
  return smod;

}


size_model controller::compute_post_model(int size, vector<int>& best_model){

  int totalc = (1<<s)-1;
  
  size_model smod;
  smod.size = size;
  smod.log10_sum_post = 0; 
 
 
  gsl_combination * gc;
  
  gc = gsl_combination_calloc (cand_set.size(), size);
  




  vector<vector<int> > mc_vec;

  


  // generate all combinations
  while(1){ 
    vector<int> cm;
    for(int i=0;i<size;i++){
      int index = cand_set[gsl_combination_get(gc,i)];
      cm.push_back(index);
    }
    
    mc_vec.push_back(cm);
    if(gsl_combination_next (gc) != GSL_SUCCESS)
      break;

  }
  gsl_combination_free(gc);
  
  
  vector<double> post_vec(mc_vec.size());
  

  int ms = mc_vec.size();

  vector<vector<vector<double> > > geno_vec = pars.geno_vec;
  vector<vector<vector<double> > > covar_vec = pars.covar_vec;
  vector<vector<double> > pheno_vec = pars.pheno_vec;
  vector<double> my_phi2_vec = phi2_vec;
  vector<double> my_omg2_vec = omg2_vec;
  double my_sslr_sigma_option = sslr_sigma_option;
  

   #pragma omp parallel for num_threads(nthread)
  for(int i=0;i<ms;i++){
    
    // empty sets to start
    vector<vector<int> > mcfg = null_config;
    map<int,int> mcfg_map = null_cfg_map; 
    vector<int> cm = mc_vec[i];
    
    string name ="";
    for(int j=0;j<size;j++){
      int index = cm[j];
      name += pars.geno_map[index];
      if(j!=size-1){
	name += "&";
      }
      
      mcfg[index]=get_config(totalc);
      mcfg_map[index] = totalc;  
     
    }

   
    SSLR local_sslr ;
    local_sslr.init(pheno_vec,geno_vec,covar_vec);
    local_sslr.set_sigma_option(0.5);
    local_sslr.set_effect_vec(my_phi2_vec,my_omg2_vec);

    double log10_post = local_sslr.compute_log10_ABF(mcfg);
    log10_post +=  compute_log10_prior(mcfg_map);
    
    #pragma omp critical 
    {
    post_vec[i] = log10_post;
    smod.post_map[name] = log10_post;
    }
  }


  double max_log10_post = -9999;
  for(int i=0;i<mc_vec.size();i++){
                     
    if(post_vec[i]>max_log10_post){
      max_log10_post = post_vec[i];
      best_model = mc_vec[i];
    }
  }
    
  vector<double> wv(post_vec.size(),1.0);
  smod.log10_sum_post = sslr.log10_weighted_sum(post_vec,wv);
  

  /*
  printf("========= %2d    %.3f  [", size, max_log10_post);
  for(int i=0;i<best_model.size();i++){
    printf("%d::",best_model[i]);
  }
  printf("]\n");
  */





  return smod;

}



void controller::backward_check(vector<int> & bm){

  if(bm.size()<2)
    return;

  int totalc = (1<<s)-1;
  

  for(int sz=2;sz<bm.size();sz++){
    gsl_combination * gc;
    gc = gsl_combination_calloc (bm.size(), sz);
    vector<vector<int> > mc_vec;
    
    // generate all combinations
    while(1){
      vector<int> cm;
      for(int i=0;i<sz;i++){
	int index = bm[gsl_combination_get(gc,i)];
	cm.push_back(index);
      }

      mc_vec.push_back(cm);
      if(gsl_combination_next (gc) != GSL_SUCCESS)
	break;
      
    }
    gsl_combination_free(gc);


    // check
    for(int i=0;i<mc_vec.size();i++){

      vector<int> cm = mc_vec[i];
      vector<vector<int> > mcfg = null_config;
      map<int,int> mcfg_map = null_cfg_map;

      // find model name
      string name ="";
      for(int j=0;j<sz;j++){
	int index = cm[j];

	mcfg[index]=get_config(totalc);
	mcfg_map[index] = totalc;

	name += pars.geno_map[index];
	if(j!=sz-1){
	  name += "&";
	}
	
      }
      
      //printf("==== %s   %.3f\n", name.c_str(), szm_vec[sz-1].post_map[name]); 
      if(szm_vec[sz-1].post_map.find(name) == szm_vec[sz-1].post_map.end()){
	szm_vec[sz-1].post_map[name] =  sslr.compute_log10_ABF(mcfg) + compute_log10_prior(mcfg_map);
      }
      
      szm_vec[sz-1].update();
    }
    
  
  }

}



void controller::summarize_approx_posterior(){
  
  vector<double>  log10_pmass_vec;
  log10_pmass_vec.push_back(compute_log10_prior(null_cfg_map));

  for(int i=0;i<szm_vec.size();i++){
    log10_pmass_vec.push_back(szm_vec[i].log10_sum_post);
  }
  double val = log10_pmass_vec[log10_pmass_vec.size()-1];

  double sum = 0;
  double ratio = 1;

  if(max_size == p){
    for(int k=log10_pmass_vec.size()+1;k<=p;k++){
      ratio *= (p-k+1)*prior_ratio/k;
      sum += ratio;
    }
    
    log10_pmass_vec.push_back(val+log10(sum));
  }

  vector<double> wv(log10_pmass_vec.size(),1.0);
  log10_pnorm  = sslr.log10_weighted_sum(log10_pmass_vec,wv);
  

  Nmodel nm;
  nm.id = "NULL";
  nm.prob = pow(10, log10_pmass_vec[0]-log10_pnorm);
  nm.size = 0;
  nm.post_score = log10_pmass_vec[0];
  nmodel_vec.push_back(nm);
  
  for(int i=0;i<szm_vec.size();i++){
    map<string, double>::iterator iter;
    for(iter = szm_vec[i].post_map.begin(); iter != szm_vec[i].post_map.end(); iter++){
      Nmodel nm;
      nm.id = iter->first;
      nm.post_score = iter->second;
      nm.prob = pow(10, iter->second-log10_pnorm);
      nm.size = i+1;
      nmodel_vec.push_back(nm);
      parse_nmodel(nm);
    }
  }
  
  std::sort(nmodel_vec.begin(),nmodel_vec.end(),sort_nmodel_dec);
  double cump = 0;
  
  double msize_mean = 0;
  double msize_var = 0;

  for(int i=0;i<nmodel_vec.size();i++){
    
    string name = nmodel_vec[i].id;
      
    size_t start_pos = 0;
    while((start_pos = name.find("&", start_pos)) != std::string::npos) {
      name.replace(start_pos, 1 , "] [");
      start_pos += 3; 
    }
    
    if( nmodel_vec[i].prob >= 1e-5){
      fprintf(outfd, "%5d   %7.4e    %d    %7.3f   [%s]\n",i+1, nmodel_vec[i].prob, nmodel_vec[i].size, nmodel_vec[i].post_score, name.c_str());
    }
    cump += nmodel_vec[i].prob;
    msize_mean += nmodel_vec[i].prob*nmodel_vec[i].size;
    msize_var  += nmodel_vec[i].prob*pow(nmodel_vec[i].size,2);
  
  }
  
  

  msize_var -= pow(msize_mean,2.0);
  if(msize_var < 0){
    msize_var = 0;
  }
  
  
  fprintf(outfd, "\nPosterior expected model size: %.3f (sd = %.3f)\n", msize_mean, sqrt(msize_var));
  fprintf(outfd, "LogNC = %.5f ( Log10NC = %.3f )", log10_pnorm/log10(exp(1)), log10_pnorm);
  
  fprintf(outfd,"\nPosterior inclusion probability\n\n");
  std::sort(nsnp_vec.begin(),nsnp_vec.end(),sort_nsnp_dec_by_ip);
  
  for(int i=0;i<nsnp_vec.size();i++){
    if(nsnp_vec[i].incl_prob<1e-3&&!output_all)
      break;
    fprintf(outfd,"%5d %10s   %8.5e    %7.3f\n",i+1, nsnp_vec[i].name.c_str(), nsnp_vec[i].incl_prob, single_log10_abfv[nsnp_vec[i].name]);
  }

  

}



void controller::parse_nmodel(Nmodel nmod){
  
  istringstream iss(nmod.id);
  string token;
  while (getline(iss, token, '&')){
    string snp_id = token;
    int index;
    if(nsnp_map.find(snp_id)!=nsnp_map.end()){
      index = nsnp_map[snp_id];
    }else{
      NSNP ns;
      ns.name = snp_id;
      ns.incl_prob = 0;
      nsnp_vec.push_back(ns);
      index = nsnp_vec.size()-1;
      nsnp_map[snp_id] = index;
    }
    nsnp_vec[index].incl_prob += nmod.prob;
    
  }
}







vector<int> controller::get_config(int c){

  vector<int> cfg;

  for(int i=0;i<s;i++){
    int r = c%2;
    cfg.push_back(r);
    c = c/2;
  }

  return cfg;
}



double *controller::get_weights(vector<double>& vec){

  double max = vec[0];
  for(int i=0;i<vec.size();i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(int i=0;i<vec.size();i++){
    sum += pow(10, (vec[i]-max));
  }
  double *pp = new double[vec.size()];
  for(int i=0;i<vec.size();i++){
    pp[i] = pow(10, (vec[i]-max))/sum;
  }
  return pp;
  
}





void size_model::update(){
  
  map<string, double>::iterator iter;
  
  vector<double> postv;
  for(iter = post_map.begin(); iter != post_map.end(); iter++){
    postv.push_back(iter->second);
  }
  
  vector<double> wv(postv.size(),1.0);
  log10_sum_post = log10_weighted_sum(postv,wv);
  
}



double log10_weighted_sum(vector<double> &vec, vector<double> &wts){


  double max = vec[0];
  for(size_t i=0;i<vec.size();i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(size_t i=0;i<vec.size();i++){
    sum += wts[i]*pow(10, (vec[i]-max));
  }

  return (max+log10(sum));
}









bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs){
  return lhs.incl_prob > rhs.incl_prob;
}




bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs){
  return lhs.prob > rhs.prob;
}


