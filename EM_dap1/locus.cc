#include "classdef.h"
#include <math.h>


  


vector<double> & Locus::EM_update(){
  
  // update prior
  
  locus_pi0 = 1;
  for(int i=0;i<snpVec.size();i++){
    double wts = alpha0[0] + alpha1[0]*snpVec[i].feature;
    snpVec[i].prior_weight = exp(wts);
    locus_pi0 *= 1.0/exp(wts);
  }


  // compute log10_lik
  vector<double> BF_vec;
  vector<double> prior_vec;
  BF_vec.push_back(0);
  prior_vec.push_back(locus_pi0);
  
  for(int i=0;i<snpVec.size();i++){ 
    BF_vec.push_back(snpVec[i].log10_BF);
    prior_vec.push_back(snpVec[i].prior_wieght*locus_pi0);
  }
  
  log10_lik = log10_weighted_sum(BF_vec, prior_vec);
  
  vector<double> loc_vec(0.0,4);
   
  for(int i=0;i<snpVec.size();i++){
    
    double val = log10(snpVec[i].prior_wieght)+ log10(locus_pi0) + snpVec[i].log10_BF - log10_lik;
    double pip = pow(10, val);
    
    if(snpVec[i].feature == 0){
      loc_vec[0] = 1-pip;
      loc_vec[1] = pip;
    }else{
      loc_vec[2] = 1-pip;
      loc_vec[3] = pip;
    }
    
  }

  return loc_vec;
}




