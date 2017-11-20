using namespace std;

#include "SSLR.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>

void SSLR::init(vector<vector<double> > & Y_in, vector<vector<vector<double> > >& Xg_in, vector<vector<vector<double> > > & Xc_in){
  
  s = Y_in.size();
  p = Xg_in[0].size();
  for(int i=0;i<s;i++){
    nv.push_back(Y_in[i].size());
    qv.push_back(Xc_in[i].size());
  }
  
    
  for(int i=0;i<s;i++){
    
    // for the intercept 
    qv[i]++;    
    gsl_matrix *Y  = gsl_matrix_calloc(nv[i],1);
    gsl_matrix *Xg = gsl_matrix_calloc(nv[i],p);
    gsl_matrix *Xc = gsl_matrix_calloc(nv[i],qv[i]);
  
    for(int j=0;j<nv[i];j++){
      gsl_matrix_set(Y,j,0,Y_in[i][j]);
    }
    
    for(int k=0;k<p;k++){
      for(int j=0;j<nv[i];j++){
	gsl_matrix_set(Xg,j,k,Xg_in[i][k][j]);
      }
    }

    if(qv[i]>1){
      for(int k=1;k<qv[i];k++){
	for(int j=0;j<nv[i];j++){
	  gsl_matrix_set(Xc,j,k,Xc_in[i][k-1][j]);
	}
      }
    }
  
    for(int j=0;j<nv[i];j++){
      gsl_matrix_set(Xc,j,0,1.0);
    }
    


    Yv.push_back(Y);
    Xgv.push_back(Xg);
    Xcv.push_back(Xc);
    mv.push_back(0);
    hv.push_back(1e-4);
  }

  
  //T = eVb = eVg_inv = Gamma = Wg =0;
  eVb = eVg_inv = Gamma = Wg = Sigma_inv=0;
  compute_common(); 
 
  compute_sigma_null();
  
}


// ================ Setting parmaeters/options ======================== //


void SSLR::set_effect_vec(vector<double> &phi2_vec_in , vector<double> &omg2_vec_in){ 

  phi2_vec = phi2_vec_in;
  omg2_vec = omg2_vec_in;
  
}


void SSLR::set_IG_prior(vector<double> & hv_in, vector<double> & mv_in){
  
  hv = hv_in;
  mv = mv_in;
}


// =============== Core Computation ================================ //


void SSLR::compute_common(){
  
  
  for(int i=0;i<s;i++){
  
    gsl_matrix *XctXc = gsl_matrix_calloc(qv[i],qv[i]);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xcv[i],Xcv[i],0,XctXc);
    gsl_matrix *XctXc_inv = gsl_matrix_calloc(qv[i],qv[i]);
 


    if(qv[i]==1)
      gsl_matrix_set(XctXc_inv,0,0,1.0/gsl_matrix_get(XctXc,0,0));
    else{
      
      
      gsl_matrix *U = gsl_matrix_calloc(qv[i],qv[i]);
      gsl_matrix *V = gsl_matrix_calloc(qv[i],qv[i]);
      gsl_vector *S = gsl_vector_calloc(qv[i]);
      gsl_vector *work = gsl_vector_calloc(qv[i]);
      
      gsl_matrix_memcpy(U,XctXc);
      
      gsl_linalg_SV_decomp(U,V,S,work);
      
      gsl_matrix *SI = gsl_matrix_calloc(qv[i],qv[i]);
      for(int j=0;j<qv[i];j++){
	double val = gsl_vector_get(S,j);
	if(val>1e-8)
	  val = 1.0/val;	
	gsl_matrix_set(SI,j,j,val);
      }
      gsl_matrix *t1 = gsl_matrix_calloc(qv[i],qv[i]);
      
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,U,SI,0,t1);
      gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t1,V,0,XctXc_inv);
      
      gsl_matrix_free(t1);
      gsl_matrix_free(SI);
      gsl_matrix_free(U);
      gsl_matrix_free(V);
      gsl_vector_free(S);
      gsl_vector_free(work);

      
      

      /*
      gsl_matrix *t = gsl_matrix_calloc(qv[i],qv[i]);
      gsl_matrix_memcpy(t,XctXc);

      int ss;
      gsl_permutation * pp = gsl_permutation_alloc(qv[i]);
      gsl_linalg_LU_decomp (t, pp, &ss);
      gsl_linalg_LU_invert (t, pp, XctXc_inv);
      gsl_permutation_free(pp);
      gsl_matrix_free(t);
      

      print_matrix(XctXc_inv,qv[i],qv[i]);
      */
      
    }
  
   
     
    gsl_matrix *t1 = gsl_matrix_calloc(qv[i],nv[i]);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XctXc_inv,Xcv[i],0,t1);
  
    gsl_matrix *t2 = gsl_matrix_calloc(nv[i],nv[i]);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Xcv[i],t1,0,t2);


  
    gsl_matrix *T = gsl_matrix_calloc(nv[i],nv[i]);
    for(int k=0;k<nv[i];k++){
      for(int j=0;j<nv[i];j++){
	double v = -gsl_matrix_get(t2,k,j);
	if(k==j)
	  v += 1;
	gsl_matrix_set(T,k,j,v);
      }
    }
   
    Tv.push_back(T);
    
    

    gsl_matrix_free(XctXc);
    gsl_matrix_free(XctXc_inv);
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
  }
}



void SSLR::compute_sigma(vector<vector<int> >& config){

  // use the null Sigma estimated initially
  
  if(sigma_vec.size()!=0)
    sigma_vec.clear();

  if(sigma_option < 1e-6){    
    sigma_vec = sigma0_vec; 
  }else{
    compute_sigma_mle(config);
  }
  
  invert_Sigma();
  return;
  
}


void SSLR::compute_sigma_null(){
  
  // compute S0, the under the null
  for(int i=0;i<s;i++){
    
    gsl_matrix *t1 = gsl_matrix_calloc(nv[i],1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Tv[i],Yv[i],0,t1);
  
    gsl_matrix *t2 = gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Yv[i],t1,0,t2);
  


    double s0 = gsl_matrix_get(t2,0,0)/nv[i];
    sigma0_vec.push_back(s0);
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
  
  }
    
  
}





void SSLR::compute_sigma_mle(vector<vector<int> >& config){
  

  for(int i=0;i<s;i++){

    gsl_vector *cpv = gsl_vector_calloc(nv[i]);
    vector<double> fac_vec;
    
    vector<int> index_vec;
    for(int j=0;j<p;j++){
      if(config[j][i] == 1){
	index_vec.push_back(j);
      }
    }
   

    int msize = qv[i]+index_vec.size();

    
    gsl_matrix *yv = gsl_matrix_calloc(nv[i],1);
    gsl_matrix_get_col(cpv,Yv[i],0);
    gsl_matrix_set_col(yv,0,cpv);
    
    gsl_matrix *Xv = gsl_matrix_calloc(nv[i],msize);
    for(int j=0;j<qv[i];j++){
      gsl_matrix_get_col(cpv,Xcv[i],j);
      gsl_matrix_set_col(Xv,j,cpv);
    }
    
    for(int j=0;j<index_vec.size();j++){
      gsl_matrix_get_col(cpv,Xgv[i],index_vec[j]);
      gsl_matrix_set_col(Xv,qv[i]+j,cpv);
    }

       
    double factor = 1;
    gsl_matrix *resv = compute_residual(yv,Xv,msize,nv[i],factor);
    //fac_vec.push_back(sqrt(factor));
    fac_vec.push_back(1);
    

    

    gsl_matrix *S =gsl_matrix_calloc(1,1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,resv,resv,0,S);
    double val = gsl_matrix_get(S,0,0)/double(mv[i]+nv[i]);
    val += hv[i]*mv[i]/(mv[i]+nv[i]);
      
    val = sigma_option*val + (1-sigma_option)*sigma0_vec[i];
    sigma_vec.push_back(val);
    
    gsl_matrix_free(yv);
    gsl_matrix_free(Xv);
    gsl_matrix_free(resv);
    gsl_vector_free(cpv);
    gsl_matrix_free(S);
    
  }
  
}
    

void SSLR::invert_Sigma(){
  
  if(Sigma_inv !=0){
    gsl_matrix_free(Sigma_inv);
  }
  
  
  Sigma_inv = gsl_matrix_calloc(s,s);
  
  for(int i=0;i<s;i++){
    gsl_matrix_set(Sigma_inv, i,i, 1.0/sigma_vec[i]);
  }
  
}




gsl_matrix *SSLR::compute_residual(gsl_matrix *y, gsl_matrix *X, int size,int n,  double &factor){
  
  gsl_matrix *XtX = gsl_matrix_calloc(size, size);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);
  
  // compute inverse of XtX (generalized inverse version)
  gsl_matrix *V = gsl_matrix_calloc(size,size);
  gsl_vector *S = gsl_vector_calloc(size);
  gsl_vector *work = gsl_vector_calloc(size);
  gsl_linalg_SV_decomp (XtX, V, S,work);
  
  gsl_matrix *t1 = gsl_matrix_calloc(size,size);
  for(int i=0;i<size;i++){
    double v = gsl_vector_get(S,i);
    if(v>1e-8){
      gsl_matrix_set(t1,i,i,1/v);
    }
  }
  gsl_matrix *t2 = gsl_matrix_calloc(size,size);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);
  
  gsl_matrix *XtX_inv = gsl_matrix_calloc(size,size);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);
  
  
  // (X'X)^{-1)X'
  gsl_matrix *t3 = gsl_matrix_calloc(size,n);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);
  
  gsl_matrix *hB = gsl_matrix_calloc(size,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,y,0,hB);

  gsl_matrix *fy = gsl_matrix_calloc(n,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,hB,0,fy);
  
  gsl_matrix *res= gsl_matrix_calloc(n,1);
  gsl_matrix_memcpy(res,y);
  gsl_matrix_sub(res,fy);
  
    

  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  
  gsl_matrix_free(XtX);
  gsl_matrix_free(XtX_inv);
  
  gsl_matrix_free(hB);
  gsl_matrix_free(fy);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);

  return res;
  
}


// construct priors according to skeleton information

void SSLR::construct_Gamma(vector<vector<int> >& indicator, vector<int> &noz_vec){
    
  if(Gamma != 0){
    gsl_matrix_free(Gamma);
    Gamma = 0;
  }
  
  Gamma = gsl_matrix_calloc(ep*s,ep*s);
  
  construct_meta_Gamma(indicator,noz_vec);
  
}




void SSLR::construct_meta_Gamma(vector<vector<int> >& indicator, vector<int> &noz_vec){
    
  for(int i=0;i<ep;i++){ 
    int index = noz_vec[i];
    gsl_matrix *t1 = gsl_matrix_calloc(s,1);
    gsl_matrix *t2 = gsl_matrix_calloc(s,s);
    
    for(int j=0;j<s;j++){
      gsl_matrix_set(t1,j,0,indicator[index][j]/sqrt(gsl_matrix_get(Sigma_inv,j,j)));    
    }
    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t1,t1,0,t2);
    for(int j=0;j<s;j++){
      for(int k=0;k<s;k++){
	gsl_matrix_set(Gamma,i*s+j,i*s+k,gsl_matrix_get(t2,j,k));
      }
    }
    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
  }
  return;
}


void SSLR::set_Wg(double phi2, double omg2){
  
  if(Wg != 0){
    gsl_matrix_free(Wg);
    Wg = 0;
  }
  
  Wg = gsl_matrix_calloc(ep*s,ep*s);  
  gsl_matrix_memcpy(Wg, Gamma);
  gsl_matrix_scale(Wg, omg2);
  
  for(int j=0;j<ep*s;j++){
    gsl_matrix_set(Wg,j,j,(omg2+phi2)*gsl_matrix_get(Gamma,j,j));
  }



  //if(omg2 > 0 ){
  gsl_matrix *P = gsl_matrix_calloc(ep*s, ep*s);

  // construct a permutation matrix
  // for (group, marker) = (i,j) 
  // prior order is given by j*s + i
  // regression order is given by p*i+j
  // need to permute (j*s+i) -> p*i+j
    
  for(int i=0;i<s;i++){
    for(int j=0;j<ep;j++){
      int index = j*s+i;
      int order = ep*i+j;
      gsl_vector *pv = gsl_vector_calloc(ep*s);
      gsl_vector_set(pv,order,1.0);
      gsl_matrix_set_row(P,index,pv);
      gsl_vector_free(pv);
    }
  }
  
  
  //print_matrix(Wg,ep*s,ep*s);
  
  gsl_matrix *t1 = gsl_matrix_calloc(ep*s, ep*s);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,P,Wg,0,t1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t1,P,0,Wg);
  
  //print_matrix(Wg,ep*s,ep*s);


  gsl_matrix_free(t1);
  gsl_matrix_free(P);
  //}
  
}


void SSLR::compute_stats(vector<int>& noz_vec){ 
  

  eVg_inv = gsl_matrix_calloc(ep*s,ep*s);
  eVb = gsl_matrix_calloc(ep*s,1);
  
  for(int i=0;i<s; i++){
    
    gsl_matrix *eXg = gsl_matrix_calloc(nv[i],ep);
    gsl_vector *mv = gsl_vector_calloc(nv[i]);
    for(int j=0;j<ep;j++){
      int index = noz_vec[j];
      gsl_matrix_get_col(mv,Xgv[i],index);
      gsl_matrix_set_col(eXg,j,mv);
    }
    
    
    
     
    gsl_matrix *G = gsl_matrix_calloc(nv[i],ep);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Tv[i],eXg,0,G);


    gsl_matrix *eKg = gsl_matrix_calloc(ep,ep);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,G,G,0,eKg);
    
    
    int offset = i*ep;
    double s2inv = gsl_matrix_get(Sigma_inv,i,i);
    
    
  
    gsl_matrix *t1 = gsl_matrix_calloc(ep,1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,G,Yv[i],0,t1);
    for(int j=0;j<ep;j++){
      for(int k=0;k<ep;k++){
	gsl_matrix_set(eVg_inv,offset+j,offset+k,gsl_matrix_get(eKg,j,k)*s2inv);
      }
      gsl_matrix_set(eVb,offset+j,0,gsl_matrix_get(t1,j,0)*s2inv);
    }  
  
    gsl_matrix_free(eXg);
    gsl_vector_free(mv);
    gsl_matrix_free(G);
    gsl_matrix_free(eKg);
    gsl_matrix_free(t1);
  }


}





double SSLR::compute_log10_ABF(gsl_matrix *Wg){
  // 1. I+Vg^-1Wg
  gsl_matrix *t1 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,eVg_inv,Wg,0,t1);
  
  for(int i=0;i<ep*s;i++){
    gsl_matrix_set(t1,i,i,gsl_matrix_get(t1,i,i)+1);
  }
 
  
  // 2.  (I+Vg^-1Wg)^-1  and determinant
  int ss;
  gsl_permutation *pp = gsl_permutation_alloc(ep*s);
  gsl_linalg_LU_decomp (t1, pp, &ss);
  double log_detVal = gsl_linalg_LU_lndet(t1);  
  gsl_matrix *t2 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_linalg_LU_invert (t1, pp, t2);
  
  
  //3. Wg(I+Vg^-1Wg)^-1
  gsl_matrix *t3 = gsl_matrix_calloc(ep*s,ep*s);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,Wg,t2,0,t3);
  
  //4. quadratic form
  gsl_matrix *t4 = gsl_matrix_calloc(1,ep*s);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,eVb,t3,0,t4);
  gsl_matrix *t5 = gsl_matrix_calloc(1,1);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t4,eVb,0,t5);
  
  double rst = .5*gsl_matrix_get(t5,0,0);


  rst += -0.5*log_detVal;
  //printf("log_detVal = %f\n",log_detVal*(-0.5));

  gsl_permutation_free(pp);

  gsl_matrix_free(t1);
  gsl_matrix_free(t2);
  gsl_matrix_free(t3);
  gsl_matrix_free(t4);
  gsl_matrix_free(t5);

  return rst/log(10.0);

}


  




double SSLR::compute_log10_ABF(vector<vector<int> >& indicator){
  
  
  // construct the skeleton from the indicator
  vector<int> noz_vec; // non_zero prior SNP index
  
  for(int i=0;i<p;i++){
    int flag = 0;
    for(int j=0;j<s;j++){
      if(indicator[i][j]!=0){
	flag=1;
	break;
      }
    } 
    if(flag==1)
      noz_vec.push_back(i);
  }
  
  ep = noz_vec.size();
  if(ep == 0)
    return 0.0;

 
  // 1. first compute Sigma
  compute_sigma(indicator);
  
 



  // 2. pre-compute necessary pieces common to all (phi2,omg2)
  compute_stats(noz_vec);
  

  // 2. scaled skeleton cov
  construct_Gamma(indicator,noz_vec); 
  
  vector<double> rst_vec;
  vector<double> wts_vec;
  
  int size = omg2_vec.size();
  
  for(int i=0;i<size;i++){

    set_Wg(phi2_vec[i], omg2_vec[i]);
    rst_vec.push_back(compute_log10_ABF(Wg));
    wts_vec.push_back(1.0/double(size));
    gsl_matrix_free(Wg);
    Wg = 0;
  }
  
  gsl_matrix_free(eVb);
  gsl_matrix_free(eVg_inv);
  gsl_matrix_free(Gamma);
  eVb=eVg_inv=Gamma=0;

  return log10_weighted_sum(rst_vec, wts_vec);

} 
		      


double SSLR::compute_log10_ABF(vector<vector<int> >& indicator, double phi2, double omg2){
 
  // construct the skeleton from the indicator
  vector<int> noz_vec; // non_zero prior SNP index
  
  for(int i=0;i<p;i++){
    int flag = 0;
    for(int j=0;j<s;j++){
      if(indicator[i][j]!=0){
	flag=1;
	break;
      }
    } 
    if(flag==1)
      noz_vec.push_back(i);
  }
  
  ep = noz_vec.size();
  if(ep == 0)
    return 0.0;


  // 1. first compute Sigma
  compute_sigma(indicator);
  // 2. pre-compute necessary pieces common to all (phi2,omg2)
  compute_stats(noz_vec);

  // 2. scaled skeleton cov
  construct_Gamma(indicator,noz_vec); 
   
  set_Wg(phi2, omg2);
  double rst = compute_log10_ABF(Wg);
  
  gsl_matrix_free(Wg);
  gsl_matrix_free(eVb);
  gsl_matrix_free(eVg_inv);
  gsl_matrix_free(Gamma);
  Wg=eVb=eVg_inv=Gamma=0;

  return rst;

} 






// utility

// simple linear regression to estimate SNP effect
vector<double> SSLR::get_SNP_effect(int p){

  vector<double> rst;
  for(int i=0;i<s;i++){
    
    double *yv = new double[nv[i]];
    double *gv = new double[nv[i]];
    for(int j=0;j<nv[i];j++){
      yv[j] = gsl_matrix_get(Yv[i],j,0);
      gv[j] = gsl_matrix_get(Xgv[i],j,p);
    }
    
    double bhat = 0;
    double muhat = 0;
    double vmu=0;
    double vb=0;
    double covmb = 0;
    double sumsq =0;
    if(gsl_stats_variance(gv,1,nv[i])>0)
      gsl_fit_linear (gv, 1, yv, 1, nv[i], &muhat, &bhat, &vmu, &covmb, &vb, &sumsq);
    rst.push_back(bhat);
    rst.push_back(vb);
    delete[] yv;
    delete[] gv;
  }

  return rst;


}



double SSLR::log10_weighted_sum(vector<double> &vec, vector<double> &wts){


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





void SSLR::print_matrix(gsl_matrix *M, int a, int b){
  
  for(int i=0;i<a;i++){
    for(int j=0;j<b;j++){
      //printf("%f  ",gsl_matrix_get(M,i,j));
      printf("%e  ",gsl_matrix_get(M,i,j));
    }
    printf("\n");
  }
  printf("\n");
}



SSLR::~SSLR(){
  
  for(int i=0;i<s;i++){
    gsl_matrix_free(Yv[i]);
    gsl_matrix_free(Xgv[i]);
    gsl_matrix_free(Xcv[i]);
    gsl_matrix_free(Tv[i]);
  }
    
  if(Wg!=0){
    gsl_matrix_free(Wg);
  }
      
  if(Sigma_inv != 0)
    gsl_matrix_free(Sigma_inv);
 
  if(Gamma!=0)
    gsl_matrix_free(Gamma);  
  
  if(eVg_inv!=0)
    gsl_matrix_free(eVg_inv);  
  
  if(eVb!=0)
    gsl_matrix_free(eVb);  
  
  hv.clear();
  mv.clear();
  sigma0_vec.clear();
  sigma_vec.clear();
  Yv.clear();
  Xgv.clear();
  Xcv.clear();
  Tv.clear();
  omg2_vec.clear();
  phi2_vec.clear();
  
}




