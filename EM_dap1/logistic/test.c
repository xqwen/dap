#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "logistic.h"


#define NOBS 6694234
#define NUM_CAT_VARS 2


int main (void)
{
  int i, j; 


  FILE *infile;

  char *inname = "../../glm.dat";

  gsl_vector_int *nlevs = gsl_vector_int_alloc(NUM_CAT_VARS); /* vector with max cat levels.*/

  gsl_vector_int *catmin = gsl_vector_int_alloc(NUM_CAT_VARS); /* vector with max cat*/
  gsl_vector_int *catmax = gsl_vector_int_alloc(NUM_CAT_VARS); /* vector with min cat*/

  gsl_vector_int *cumlevs = gsl_vector_int_alloc(NUM_CAT_VARS);


  gsl_matrix_int *m = gsl_matrix_int_alloc (NOBS, NUM_CAT_VARS); /* matrix with categ. vars.*/
  gsl_vector *y = gsl_vector_alloc(NOBS); /* vector of probs. to predict  */
  gsl_vector *yhat = gsl_vector_alloc(NOBS); /* vector of probs. predicted  */

  gsl_vector *params;

  double mLogLik=0;
  double mse=0;

  // Reading data
  infile = fopen(inname, "r");
  if (!infile) {
    printf("Couldn't open %s for reading\n",inname);
    return 0;
  }

  for (i = 0; i < NOBS; i++){
    fscanf(infile,"%*s%*s%lf",gsl_vector_ptr(y,i));
    for (j = 0; j < NUM_CAT_VARS; j++)
    {
      fscanf(infile,"%d",gsl_matrix_int_ptr(m, i, j));
      if(gsl_matrix_int_get(m,i,j)<gsl_vector_int_get(catmin,j))
	gsl_vector_int_set(catmin,j,gsl_matrix_int_get(m,i,j));
      if(gsl_matrix_int_get(m,i,j)>gsl_vector_int_get(catmax,j))
	gsl_vector_int_set(catmax,j,gsl_matrix_int_get(m,i,j));
    }
  }

  // Rescaling factors so they are in the 0,NumLevs range
  for (j = 0; j < NUM_CAT_VARS; j++){
    gsl_vector_int_set(nlevs,j,gsl_vector_int_get(catmax,j)-gsl_vector_int_get(catmin,j)+1);

    for (i = 0; i < NOBS; i++)
      gsl_matrix_int_set(m,i,j,
			 gsl_matrix_int_get(m,i,j)-gsl_vector_int_get(catmin,j));
  }

  gsl_vector_int_set(cumlevs,0,1+gsl_vector_int_get(nlevs,0)-1);
  for (j = 1; j < NUM_CAT_VARS; j++)
    gsl_vector_int_set(cumlevs,j,
		       gsl_vector_int_get(cumlevs,j-1)+gsl_vector_int_get(nlevs,j)-1);

  for (j = 0; j < NUM_CAT_VARS; j++){
    printf("#Factor var(%d):%d,%d,%d,%d\n"
	   ,j
	   ,gsl_vector_int_get(catmin,j)
	   ,gsl_vector_int_get(catmax,j)
	   ,gsl_vector_int_get(nlevs,j)
	   ,gsl_vector_int_get(cumlevs,j));
  }

  // Setting up input parameter vector intitialized to zeros
  params=gsl_vector_alloc(gsl_vector_int_get(cumlevs,NUM_CAT_VARS-1));  

  // Playing with the variables
  //  gsl_vector_set(params,0,-8);
  // gsl_vector_set(params,100,0);
  //gsl_vector_set(params,201,0.5);
  //gsl_vector_set(params,202,0.9);



  // Testing fLogit function with intiial parameters:
  mLogLik=fLogitCat(params,m,nlevs,y,0.0,0.0);
  printf("#Initial -log(Lik(0))=%lf\n",mLogLik);

  logistic_pred(params,m,nlevs,yhat);
  mse=0;
  for (i = 0; i < NOBS; i++){
    double aux=(gsl_vector_get(y,i)-gsl_vector_get(yhat,i));
    mse+=aux*aux;
  }
  mse=mse/(double)NOBS;
  printf("#MSE(y-yhat)_start:%lg\n",mse);


  // OPTIMIZER
  logistic_fit(params,m,nlevs,y,0.0,0.0);
  
  logistic_pred(params,m,nlevs,yhat);

  mse=0;
  for (i = 0; i < NOBS; i++){
    double aux=(gsl_vector_get(y,i)-gsl_vector_get(yhat,i));
    mse+=aux*aux;
  }
  mse=mse/(double)NOBS;
  printf("#MSE(y-yhat)_end:%lg\n",mse);

  // Print top 20 lines
  printf ("#test order:"); 
  for (j = 0; j < 40; j++)
    printf ("%d,",m->data[j]);
  printf ("\n");
  

  // Debug top 20 lines 
  for (i = 0; i < 20; i++){  
    printf("%lg",gsl_vector_get(y,i));
    for (j = 0; j < NUM_CAT_VARS; j++)
      printf ("\t%d",gsl_matrix_int_get (m, i, j));
    printf("\n");
  }
  // 
 
  gsl_matrix_int_free (m);
  gsl_vector_free (y);
  gsl_vector_int_free(nlevs);
  gsl_vector_int_free(cumlevs);
  gsl_vector_int_free(catmin);
  gsl_vector_int_free(catmax);

  gsl_vector_free(params);

  return 0;
}
