using namespace std;
#include "SSLR_SS.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>

void SSLR_SS::init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_){
  

  
  n = n_;
  p = GtG_->size1;
  yty = yty_;
  
  GtG = gsl_matrix_calloc(p,p);
  Gty = gsl_matrix_calloc(p,1);
  gsl_matrix_memcpy(GtG, GtG_);
  gsl_matrix_memcpy(Gty, Gty_);
  
}



SSLR_SS::~SSLR_SS(){
 
  gsl_matrix_free(GtG);
  gsl_matrix_free(Gty);

}

// ================ Setting parmaeters/options ======================== //


void SSLR_SS::set_effect_vec(vector<double> &phi2_vec_){ 

  phi2_vec = phi2_vec_;
  
}



double SSLR_SS::compute_log10_ABF(vector<vector<int> > & indicator){
 
  vector<int> single_cfg;
  for(int i=0;i<indicator.size();i++){
    single_cfg.push_back(indicator[i][0]);
  }
  return(compute_log10_ABF(single_cfg));

}


double SSLR_SS::compute_log10_ABF(vector<int> & indicator){
 

  // construct the sub-matrices
  
  vector<double> rstv;
  vector<double> wv;
  
   
  int ep = 0;
  int count = 0;
  std::map<int,int> imap;


  for(int i=0;i<indicator.size();i++){
    if(indicator[i]==1){
      ep++;
      imap[i] = count++;
    }
  }

  gsl_matrix *XtX = gsl_matrix_calloc(ep,ep);
  gsl_matrix *Xty = gsl_matrix_calloc(ep,1);
 

  for(int i=0;i<p;i++){
    if(indicator[i] == 0)
      continue;
    gsl_matrix_set(Xty, imap[i],0,gsl_matrix_get(Gty,i,0));
    for(int j=0;j<p;j++){
      if(indicator[j] == 1){
	
	double val = gsl_matrix_get(GtG,i,j);
	gsl_matrix_set(XtX,imap[i], imap[j],val);
      }
    }
  }

  
  // compute inverse of XtX (generalized inverse version)
  gsl_matrix *V = gsl_matrix_calloc(ep,ep);
  gsl_vector *S = gsl_vector_calloc(ep);
  gsl_vector *work = gsl_vector_calloc(ep);
  gsl_linalg_SV_decomp (XtX, V, S,work);


  
  gsl_matrix *t1 = gsl_matrix_calloc(ep,ep);
  for(int i=0;i<ep;i++){
    double v = gsl_vector_get(S,i);
    if(v>1e-8){
      gsl_matrix_set(t1,i,i,1/v);
    }
  }
  gsl_matrix *t2 = gsl_matrix_calloc(ep,ep);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);


  
  gsl_matrix *XtX_inv = gsl_matrix_calloc(ep,ep);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);


  
  
  // compute sigma
  gsl_matrix *t3 = gsl_matrix_calloc(ep,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,XtX_inv,Xty,0,t3);
  
  gsl_matrix *t4 = gsl_matrix_calloc(1,1);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xty,t3,0,t4);
  
  double res = yty - gsl_matrix_get(t4,0,0);
 
  double sigma2 = (sigma_option*res + (1-sigma_option)*yty)/n;
  
 
  // compute BF
  for(int i=0;i<phi2_vec.size();i++){

    double det = 1;
    double phi2 = phi2_vec[i];
   
    gsl_matrix *tt1 = gsl_matrix_calloc(ep,ep);
    
    for(int j=0;j<ep;j++){
      det = det*(1+phi2*gsl_vector_get(S,j));

      double v = gsl_vector_get(S,j)+1.0/phi2_vec[i];
      if(v>1e-8){
	gsl_matrix_set(tt1,j,j,1.0/v);
      }
    }
    
    
    

    gsl_matrix *tt2 = gsl_matrix_calloc(ep,ep);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,tt1,0,tt2);

    gsl_matrix *IWV = gsl_matrix_calloc(ep,ep);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,tt2,V,0,IWV);
    
    gsl_matrix *tt3 = gsl_matrix_calloc(ep,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,IWV,Xty,0,tt3);

    gsl_matrix *tt4 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xty,tt3,0,tt4);
    
    double log10BF = (-0.5*log(fabs(det))+0.5*gsl_matrix_get(tt4,0,0)/sigma2)/log(10);
    
    rstv.push_back(log10BF);
    wv.push_back(1.0/phi2_vec.size());
  
    
    gsl_matrix_free(tt1);
    gsl_matrix_free(tt2);
    gsl_matrix_free(IWV);
    gsl_matrix_free(tt3);
    gsl_matrix_free(tt4);    


  }
  
  
  
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(XtX_inv);

  gsl_matrix_free(t3);
  gsl_matrix_free(t4);

  gsl_matrix_free(XtX);
  gsl_matrix_free(Xty);


  double rst =  log10_weighted_sum(rstv,wv);
  
  return rst;
} 






// utility

double SSLR_SS::log10_weighted_sum(vector<double> &vec, vector<double> &wts){


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





void SSLR_SS::print_matrix(gsl_matrix *M, int a, int b){
  
  for(int i=0;i<a;i++){
    for(int j=0;j<b;j++){
      //printf("%f  ",gsl_matrix_get(M,i,j));
      printf("%e  ",gsl_matrix_get(M,i,j));
    }
    printf("\n");
  }
  printf("\n");
}



