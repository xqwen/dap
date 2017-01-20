using namespace std;

#include "controllerZ.h"
#include <fstream>
#include <sstream>
#include <math.h>



void controllerZ::set_subset(char *snp_list){
  
  ifstream sfile(snp_list);
  string line;
  istringstream ins;
   
  set_filter = 1;
  
 while(getline(sfile,line)){

    ins.clear();
    ins.str(line);

    string snp;
    if(ins>>snp){
      use_snp[snp] = 1;
    }
 }
 
 sfile.close();

}


void controllerZ::init(char *zval_file){
  
  ifstream gfile(zval_file);
  string line;
  istringstream ins;
  int count = 0;
  while(getline(gfile,line)){

    ins.clear();
    ins.str(line);

    string snp;
    string loc;
    double zval;
    if(ins>>snp>>loc>>zval){
      
      if(set_filter && use_snp[snp] != 1){
	continue;
      }
      
      geno_rmap[snp] = count;
      geno_map[count] = snp;
      count++;
      log10_BF_vec.push_back(compute_log10_BF(zval));
    }
  }
  
  p = int(log10_BF_vec.size());
  
  gfile.close();
  

}


void controllerZ::init_BF(char *bf_file){

  ifstream gfile(bf_file);
  string line;
  istringstream ins;
  int count = 0;
  while(getline(gfile,line)){

    ins.clear();
    ins.str(line);

    string snp;
    string loc;
    double log10_bf;
    if(ins>>snp>>loc>>log10_bf){

      if(set_filter && use_snp[snp] != 1){
        continue;
      }

      geno_rmap[snp] = count;
      geno_map[count] = snp;
      count++;
      log10_BF_vec.push_back(log10_bf);
    }
  }

  p = int(log10_BF_vec.size());
  gfile.close();
  

}





void controllerZ::set_prior(double pi1){
  

  if(pi1<=0 || pi1>=1.0){
    pi1 = 1.0/double(p);
  }
  
  for (int i=0;i<p;i++){
    pvec.push_back(pi1);
  }
  return;
}



void controllerZ::set_prior(char *prior_file){
  
  if(prior_file == 0){
    set_prior(-1.0);
    return;
  }

  
  ifstream ifile(prior_file);
  string line;
  istringstream ins;

  for(int i=0;i<p;i++){
    pvec.push_back(0.0);
  }

  while(getline(ifile,line)){

    ins.clear();
    ins.str(line);
    double value;
    string snp_name;
    ins>>snp_name;

    if(set_filter && use_snp[snp_name] != 1){
      continue;
    }



    int index = geno_rmap[snp_name];
    ins>>value;

    pvec[index] = value;
  }
  return;
}
    



void controllerZ::run(){
  
  vector<double> wts(p+1,1.0);
  vector<double> postv;
  for(int i=0;i<p;i++){
    double val = log10(pvec[i])-log10(1.0 - pvec[i])+ log10_BF_vec[i];
    postv.push_back(val);
  }

  postv.push_back(0.0);
  double log10NC = log10_weighted_sum(postv,wts);
  
  for(int i=0;i<p;i++){
    string snp = geno_map[i];
    double pip = pow(10.0,postv[i]-log10NC);
    printf("%15s    %9.3e   %7.3f\n",snp.c_str(), pip,log10_BF_vec[i]);
    
  }

}


double compute_log10_BF(double z_score){
  double kv[4] = {1,4,16,25};
  int size = 4;
  double z2 = pow(z_score, 2.0);
  vector<double> rstv;
  vector<double> wtv(size,1.0/double(size));
  for(int i=0;i<size;i++){
    
    double val = 0.5*log(1/(1+kv[i])) + 0.5*z2*(kv[i]/(1+kv[i]));
    rstv.push_back(val/log(10));
  }

  return log10_weighted_sum(rstv,wtv);
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


