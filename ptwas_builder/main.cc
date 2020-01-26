using namespace std;

#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "parser.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>

#define LENGTH 1280

void parse_fm(char *fm_file);
double run_regression(vector<int> &index_vec, double prob);
double single_snp_analysis(int index);


vector<double> prob_vec;
vector<vector<string> > model_vec;


vector<vector<double> > genov;
gsl_matrix *y;
int n;
vector<double> coef_vec;
vector<double> z2_vec;

string gene;
string tissue;
double weight_thresh;

int main (int argc, char **argv){ 



    char data_file[LENGTH];
    char fm_file[LENGTH];
    char tissue_name[LENGTH];
    char gene_name[LENGTH];
    char output_file[LENGTH];

    memset(data_file, 0, LENGTH);
    memset(fm_file,0,LENGTH);
    memset(tissue_name,0,LENGTH);
    memset(gene_name,0, LENGTH);
    memset(output_file, 0, LENGTH);

    int er2_only = 0;
    weight_thresh = 0;

    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
            strcpy(data_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-f")==0 || strcmp(argv[i], "-fm")==0){
            strcpy(fm_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "--er2")==0){
            er2_only = 1;
            continue;
        }


        if(strcmp(argv[i], "-tissue")==0 || strcmp(argv[i], "-t")==0){
            strcpy(tissue_name, argv[++i]);
            continue;
        }



        if(strcmp(argv[i], "-gene")==0 || strcmp(argv[i], "-g")==0){
            strcpy(gene_name, argv[++i]);
            continue;
        }
        
        if(strcmp(argv[i], "-thresh")==0 || strcmp(argv[i], "-wt")==0){
            weight_thresh = fabs(atof(argv[++i]));
        }

        if(strcmp(argv[i], "-o")==0 ){
            strcpy(output_file, argv[++i]);
        }

    }

    if(strlen(tissue_name)==0){
        tissue = string("TISSUE");
    }else{
        tissue = string(tissue_name);
    }

    if(strlen(gene_name)==0){
        fprintf(stderr, "Error: gene name is not specified\n");
    }else{
        gene = string(gene_name);
    }



    parse_fm(fm_file);
    // parsing data file
    parser pars;
    pars.process_data(data_file);
    
    genov = pars.geno_vec[0];
    for(int i=0;i<genov.size();i++){
        coef_vec.push_back(0);
    }


    n = pars.pheno_vec[0].size();

    y = gsl_matrix_calloc(n,1);
    for(int i=0;i<n;i++){
        gsl_matrix_set(y,i,0, pars.pheno_vec[0][i]);
    }


    //single snp analysis
    vector<double> z2_vec;
    for(int i=0; i<coef_vec.size();i++){
        z2_vec.push_back(single_snp_analysis(i));
    }


    // run regression  
    double cump = 0;
    double ER2 = 0;
    for(int i = 0; i< model_vec.size(); i++){
        vector<string> model = model_vec[i];
        vector<int> index;

        cump += prob_vec[i];

        if(model.size() == 0){
            continue;
        }

        for(int k =0 ;k<model.size();k++){
            index.push_back(pars.geno_rmap[model[k]]);
        }

        ER2  += run_regression(index, prob_vec[i]);
    }

    FILE *out_fd = stdout;
    if(strlen(output_file)>0){
        out_fd = fopen(output_file, "w");
    }
    if(!er2_only){


        for(int i=0;i<coef_vec.size();i++){
            if(fabs(coef_vec[i]) >= weight_thresh)
                fprintf(out_fd, "%s\t%s\t%s\t%7.3e\t%9.3e\t%7.3e\n",tissue.c_str(), gene.c_str(), pars.geno_map[i].c_str(), coef_vec[i],z2_vec[i],2*gsl_cdf_ugaussian_P(-fabs(z2_vec[i])));
        }

    }else{
        printf("%s  %7.3e\n",fm_file, ER2);
    }


}


void parse_fm(char *grid_file){

    ifstream gfile(grid_file);
    string line;
    istringstream ins;


    while(getline(gfile,line)){

        if(line.find("[") == string::npos){
            continue;
        }
        ins.clear();
        ins.str(line);
        int index;
        double prob;
        int size;
        double pscore;

        ins>>index>>prob>>size>>pscore;
        vector<string> model;
        prob_vec.push_back(prob);
        string ms;
        for(int i=0;i<size;i++){
            ins>>ms;
            int length = ms.length();
            string snp = ms.substr(1,length-2);
            model.push_back(snp);
        }
        model_vec.push_back(model);
    }

}





double run_regression(vector<int> &index_vec, double weight){

    int p = index_vec.size();

    gsl_matrix *X = gsl_matrix_calloc(n, p+1);
    for(int i=0;i<n;i++){
        gsl_matrix_set(X,i,0,1.0);
    }
    for(int i=0;i<p;i++){
        vector<double> vec = genov[index_vec[i]];
        for(int j=0;j<n;j++){
            gsl_matrix_set(X,j,i+1,vec[j]);
        }
    }
    /*
       for(int i=0;i<n;i++){
       printf("%f  ", gsl_matrix_get(y,i,0));
       for(int j=1;j<p+1;j++){
       printf("%f ", gsl_matrix_get(X,i,j));
       }
       printf("\n");
       }	       
       */
    gsl_matrix *XtX = gsl_matrix_calloc(p+1,p+1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);

    // compute inverse of XtX (generalized inverse version)
    gsl_matrix *V = gsl_matrix_calloc(p+1,p+1);
    gsl_vector *S = gsl_vector_calloc(p+1);
    gsl_vector *work = gsl_vector_calloc(p+1);
    gsl_linalg_SV_decomp (XtX, V, S,work);

    gsl_matrix *t1 = gsl_matrix_calloc(p+1,p+1);
    for(int i=0;i<p+1;i++){
        double v = gsl_vector_get(S,i);
        if(v>1e-8){
            gsl_matrix_set(t1,i,i,1/v);
        }
    }
    gsl_matrix *t2 = gsl_matrix_calloc(p+1,p+1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);

    gsl_matrix *XtX_inv = gsl_matrix_calloc(p+1,p+1);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);

    // (X'X)^{-1)X'
    gsl_matrix *t3 = gsl_matrix_calloc(p+1,n);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);


    gsl_matrix *bv = gsl_matrix_calloc(p+1,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,y,0,bv);
    for(int i=0;i<p;i++){
        int k = index_vec[i];
        coef_vec[k] += gsl_matrix_get(bv,i+1,0)*weight;
    }

    gsl_matrix *Yhat = gsl_matrix_calloc(n,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,bv,0,Yhat);

    double resid = 0;
    double meany = 0;
    for(int i=0;i<n;i++){
        double yh = gsl_matrix_get(Yhat,i,0);
        double yo = gsl_matrix_get(y,i,0);
        resid += pow(yh-yo,2);
        meany += yo;
    }

    meany = meany/n;

    double syy = 0;
    for(int i=0;i<n;i++){
        syy += pow(gsl_matrix_get(y,i,0) - meany,2);
    }

    double er2 = weight*(1-resid/syy);

    /*
       double s2 = resid/(n-p-1);

       for(int i=0;i<p;i++){
       int k = index_vec[i];
       var_vec[k] += s2*gsl_matrix_get(XtX_inv, i+1, i+1)*weight*weight;
       }
       */

    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
    gsl_matrix_free(t3);

    gsl_matrix_free(XtX);
    gsl_matrix_free(XtX_inv);
    gsl_matrix_free(X);

    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(bv);

    gsl_matrix_free(Yhat);
    return er2;
}



double single_snp_analysis(int index){

    gsl_matrix *X = gsl_matrix_calloc(n, 2);
    for(int i=0;i<n;i++){
        gsl_matrix_set(X,i,0,1.0);
    }
    vector<double> vec = genov[index];
    for(int j=0;j<n;j++){
        gsl_matrix_set(X,j,1,vec[j]);
    }

    gsl_matrix *XtX = gsl_matrix_calloc(2,2);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0,XtX);

    // compute inverse of XtX (generalized inverse version)
    gsl_matrix *V = gsl_matrix_calloc(2,2);
    gsl_vector *S = gsl_vector_calloc(2);
    gsl_vector *work = gsl_vector_calloc(2);
    gsl_linalg_SV_decomp (XtX, V, S,work);

    gsl_matrix *t1 = gsl_matrix_calloc(2,2);
    for(int i=0;i<2;i++){
        double v = gsl_vector_get(S,i);
        if(v>1e-8){
            gsl_matrix_set(t1,i,i,1/v);
        }
    }
    gsl_matrix *t2 = gsl_matrix_calloc(2,2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);

    gsl_matrix *XtX_inv = gsl_matrix_calloc(2,2);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,XtX_inv);

    // (X'X)^{-1)X'
    gsl_matrix *t3 = gsl_matrix_calloc(2,n);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,XtX_inv,X,0,t3);


    gsl_matrix *bv = gsl_matrix_calloc(2,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,t3,y,0,bv);
    gsl_matrix *Yhat = gsl_matrix_calloc(n,1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,X,bv,0,Yhat);

    double resid = 0;
    for(int i=0;i<n;i++){
        double yh = gsl_matrix_get(Yhat,i,0);
        double yo = gsl_matrix_get(y,i,0);
        resid += pow(yh-yo,2);
    }

    double s2 = resid/(n-2);
    double var = s2*gsl_matrix_get(XtX_inv, 1, 1);
    double z = 0;
    if(var>0){
        z = gsl_matrix_get(bv,1,0)/sqrt(var);
    }


    gsl_matrix_free(t1);
    gsl_matrix_free(t2);
    gsl_matrix_free(t3);

    gsl_matrix_free(XtX);
    gsl_matrix_free(XtX_inv);
    gsl_matrix_free(X);

    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_matrix_free(bv);

    gsl_matrix_free(Yhat);
    return z;
}
