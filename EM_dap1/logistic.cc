#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>
#include "logistic.h"

// I need to bundle all the data that goes to the function to optimze together. 
typedef struct{
  gsl_matrix_int *X;
  gsl_vector_int *nlev;
  gsl_vector *y;
  double lambdaL1;
  double lambdaL2;
}fix_parm_T;


double fLogitCat(gsl_vector *beta
		 ,gsl_matrix_int *X
		 ,gsl_vector_int *nlev
		 ,gsl_vector *y
		 ,double lambdaL1
		 ,double lambdaL2)
{
  int n = y->size; 
  //  int k = X->size2; 
  int npar = beta->size; 
  double total = 0;
  double aux = 0;
/*   omp_set_num_threads(ompthr); */
/*   /\* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*\/ */
/*   /\*#pragma omp parallel for reduction (+:total)*\/ */
  for(int i = 1; i < npar; ++i)
    total += beta->data[i]*beta->data[i];
  total = (-total*lambdaL2/2);
/*   /\*#pragma omp parallel for reduction (+:aux)*\/ */
  for(int i = 1; i < npar; ++i)
    aux += (beta->data[i]>0 ? beta->data[i] : -beta->data[i]);
  total = total-aux*lambdaL1;
/* #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y) reduction (+:total) */
  for(int i = 0; i < n; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    total += y->data[i]*Xbetai-gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
} 


void logistic_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		  ,gsl_matrix_int *X  //Matrix Nobs x K 
		  ,gsl_vector_int *nlev // Vector with number categories
		  ,gsl_vector *yhat //Vector of prob. predicted by the logistic
		  )
{
  for(int i = 0; i < X->size1; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    yhat->data[i]=1/(1 + gsl_sf_exp(-Xbetai));
  }
}


/* The gradient of f, df = (df/dx, df/dy). */
void 
wgsl_optim_df (const gsl_vector *beta, void *params, 
       gsl_vector *out)
{
  fix_parm_T *p = (fix_parm_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int npar = beta->size; 
  // Intitialize gradient out necessary?
  for(int i = 0; i < npar; ++i) 
    out->data[i]= 0; 
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0 */
  for(int i = 1; i < npar; ++i)
    out->data[i]= p->lambdaL2*beta->data[i]; 
  for(int i = 1; i < npar; ++i)
    out->data[i]+= p->lambdaL1*((beta->data[i]>0)-(beta->data[i]<0));
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm];
      iParm+=p->nlev->data[k]-1;
    }
    //    total += y->data[i]*Xbetai-log(1+gsl_sf_exp(Xbetai));
    pn= -( p->y->data[i] - 1/(1 + gsl_sf_exp(-Xbetai)) );

    out->data[0]+= pn;
    iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	out->data[gsl_matrix_int_get(p->X,i,k)-1+iParm]+=pn;
      iParm+=p->nlev->data[k]-1;
    }
  }
}


/* The Hessian of f */
void 
wgsl_optim_hessian (const gsl_vector *beta, void *params, 
       gsl_matrix *out)
{
  fix_parm_T *p = (fix_parm_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int npar = beta->size; 
  // Intitialize Hessian out necessary ???
  gsl_matrix_set_zero(out);
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for(int i = 1; i < npar; ++i)
    gsl_matrix_set(out,i,i,(p->lambdaL2));  // Double check this
// L1 penalty not working yet, as not differentiable, I may need to do coordinate descent (as in glm_net)

  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double aux=0;
    double Xbetai=beta->data[0];
    int iParm2=1;
    int iParm1=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm1];
      iParm1+=p->nlev->data[k]-1;  //-1?
    }
    //    total += y->data[i]*Xbetai-log(1+gsl_sf_exp(Xbetai));
    pn= 1/(1 + gsl_sf_exp(-Xbetai));
// Add a protection for pn very close to 0 or 1?
    aux=pn*(1-pn);
    *gsl_matrix_ptr(out,0,0)+=aux;
    iParm2=1;
    for(int k2 = 0; k2 < K; ++k2) {
      if(gsl_matrix_int_get(p->X,i,k2)>0)
	*gsl_matrix_ptr(out,0,gsl_matrix_int_get(p->X,i,k2)-1+iParm2)+=aux;
      iParm2+=p->nlev->data[k2]-1;   //-1?
    }
    iParm1=1;
    for(int k1 = 0; k1 < K; ++k1) {
      if(gsl_matrix_int_get(p->X,i,k1)>0)
	*gsl_matrix_ptr(out,gsl_matrix_int_get(p->X,i,k1)-1+iParm1,0)+=aux;
      iParm2=1;
      for(int k2 = 0; k2 < K; ++k2) {
	if((gsl_matrix_int_get(p->X,i,k1)>0) && (gsl_matrix_int_get(p->X,i,k2)>0))
	  *gsl_matrix_ptr(out
			  ,gsl_matrix_int_get(p->X,i,k1)-1+iParm1
			  ,gsl_matrix_int_get(p->X,i,k2)-1+iParm2
			  )+=aux;
	iParm2+=p->nlev->data[k2]-1;  //-1?
      }
      iParm1+=p->nlev->data[k1]-1; //-1?
    }
  }
}


double wgsl_optim_f(gsl_vector *v, void *params)
{
  double mLogLik=0;
  fix_parm_T *p = (fix_parm_T *)params;
  mLogLik = fLogitCat(v,p->X,p->nlev,p->y,p->lambdaL1,p->lambdaL2);
  return mLogLik; 
}


/* Compute both f and df together. */
void 
wgsl_optim_fdf (gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{
  *f = wgsl_optim_f(x, params); 
  wgsl_optim_df(x, params, df);
}


int logistic_fit(gsl_vector *beta
		 ,gsl_matrix_int *X
		 ,gsl_vector_int *nlev
		 ,gsl_vector *y
		 ,double lambdaL1
		 ,double lambdaL2)
{

  double mLogLik=0;
  fix_parm_T p;
  int npar = beta->size; 
  int iter=0;
  double maxchange=0;

  //Intializing fix parameters
  p.X=X;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;
  
  //Initial fit
  //#ifdef _RPR_DEBUG_
  mLogLik = wgsl_optim_f(beta,&p);
  //fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
  //#endif //_RPR_DEBUG

  gsl_matrix *myH = gsl_matrix_alloc(npar,npar); /* Hessian*/
  gsl_vector *stBeta = gsl_vector_alloc(npar); /* Move */

  gsl_vector *myG = gsl_vector_alloc(npar); /* Gradient*/
  gsl_vector *tau = gsl_vector_alloc(npar); /* tau for QR*/

  for(iter=0;iter<100;iter++){ 
    wgsl_optim_hessian(beta,&p,myH);
    wgsl_optim_df(beta,&p,myG);
    gsl_linalg_QR_decomp(myH,tau);
    gsl_linalg_QR_solve(myH,tau,myG,stBeta);
    gsl_vector_sub(beta,stBeta);
    
    maxchange=0;
    for(int i=0;i<npar; i++)
      if(maxchange<fabs(stBeta->data[i]))
	maxchange=fabs(stBeta->data[i]);

#ifdef _RPR_DEBUG_
      mLogLik = wgsl_optim_f(beta,&p);
      fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,maxchange);
#endif //_RPR_DEBUG

    if(maxchange<1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
    fprintf(stderr,"#par_%d= %lf\n",i,beta->data[i]);
#endif //_RPR_DEBUG

  //Final fit
  mLogLik = wgsl_optim_f(beta,&p);
  //fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange %lf\n",iter,mLogLik,maxchange);

  gsl_vector_free (tau);
  gsl_vector_free (stBeta);
  gsl_vector_free (myG);
  gsl_matrix_free (myH);

  return 0;
}



/*


int logistic_fit_old(gsl_vector *beta
		 ,gsl_matrix_int *X
		 ,gsl_vector_int *nlev
		 ,gsl_vector *y
		 ,double lambdaL1
		 ,double lambdaL2)
{
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  double mLogLik=0;
  fix_parm_T p;

  //  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  //Intializing fix parameters
  p.X=X;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;

  //gsl minimizer set up
  my_func.n = beta->size;
  my_func.f = wgsl_optim_f;
  my_func.df = wgsl_optim_df;
  my_func.fdf = wgsl_optim_fdf;
  my_func.params = &p;
  
  //Initial fit
  mLogLik = wgsl_optim_f(beta,&p);
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);

  //Minimizer type
  //  T = gsl_multimin_fdfminimizer_conjugate_fr;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, beta->size);

  gsl_multimin_fdfminimizer_set (s, &my_func, beta, 0.5, 0.1);

  fprintf(stderr,"#s is a '%s' minimizer\n", 
	  gsl_multimin_fdfminimizer_name (s));

  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;
    
    status = gsl_multimin_test_gradient (s->gradient, 1e-5);

    mLogLik = wgsl_optim_f(s->x,&p);
    fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,s->f);
    
    if (status == GSL_SUCCESS){
      fprintf(stderr,"#Minimum found after %d iterations:\n",(int)iter);   
      for (int i = 0; i < s->x->size; i++)
	fprintf(stderr,"#par_%d= %lf\n",i,s->x->data[i]);
    }
  }while (status == GSL_CONTINUE && iter < 1000);

  for (int i = 0; i < s->x->size; i++)
    fprintf(stderr,"#par_%d= %lf\n",i,s->x->data[i]);

  for (int i = 0; i < s->x->size; i++)
    beta->data[i]=s->x->data[i];

  //Final fit
  mLogLik = wgsl_optim_f(beta,&p);
  fprintf(stderr,"#Final -log(Lik(0))=%lf\n",mLogLik);

  gsl_multimin_fdfminimizer_free (s);
  //  gsl_vector_free (x);

  return 0;
}
*/

