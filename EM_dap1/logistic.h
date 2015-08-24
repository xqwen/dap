#ifndef LOGISTIC_H_   /* Include guard */
#define LOGISTIC_H_

void logistic_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		  ,gsl_matrix_int *X  //Matrix Nobs x K 
		  ,gsl_vector_int *nlev // Vector with number categories
		  ,gsl_vector *yhat //Vector of prob. predicted by the logistic
		   );
		 
int logistic_fit(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		 ,gsl_matrix_int *X  //Matrix Nobs x K 
		 ,gsl_vector_int *nlev // Vector with number categories
		 ,gsl_vector *y //Vector of prob. to predict
		 ,double lambdaL1 // Regularization L1 0.0 if not used
		 ,double lambdaL2); // Regularization L2 0.0 if not used


double fLogitCat(gsl_vector *beta
		 ,gsl_matrix_int *X
		 ,gsl_vector_int *nlev
		 ,gsl_vector *y
		 ,double lambdaL1
		 ,double lambdaL2);

#endif // LOGISTIC_H_
