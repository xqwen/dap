using namespace std;
#include "MLR.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>

void MLR::init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_){

    n = n_;
    p = GtG_->size1;
    yty = yty_;

    GtG = gsl_matrix_calloc(p,p);
    Gty = gsl_matrix_calloc(p,1);
    gsl_matrix_memcpy(GtG, GtG_);
    gsl_matrix_memcpy(Gty, Gty_);
    use_ss = 0;
}

void MLR::init_ss(const gsl_matrix *Z_, const gsl_matrix *R_){

    p = R_->size1;
    use_ss = 1;

    R = gsl_matrix_calloc(p,p);
    Z = gsl_matrix_calloc(p,1);
    gsl_matrix_memcpy(R, R_);
    gsl_matrix_memcpy(Z, Z_);

}




void MLR::init(vector<double> &Y_in, vector<vector<double> >  &G_in){

    n = Y_in.size();
    p = G_in.size();

    yty = 0;
    Y = gsl_matrix_calloc(n,1);
    for(int i=0;i<n;i++){
        double val = Y_in[i];
        yty += pow(val,2);
        gsl_matrix_set(Y,i,0,val);
    }
    
    G = gsl_matrix_calloc(n,p);

    for(int j=0;j<p;j++){
        for(int i=0;i<n;i++){
            gsl_matrix_set(G,i,j,G_in[j][i]);
        }   
    }

    // On 2/21/2023, for computational efficiency, we no longer pre-compute G'G for individual-level data
   
   GtG =0;
   Gty = 0;
   use_ss = 0;
}



void MLR::copy(const MLR& mlr){
    if(mlr.use_ss == 0)
        copy_full(mlr);
    else
        copy_ss(mlr);
}



void MLR::copy(const MLR& mlr, vector<int> &indicator){
    if(mlr.use_ss == 0)
        copy_full(mlr, indicator);
    else
        copy_ss(mlr, indicator);
}




void MLR::copy_ss(const MLR & mlr){

    if(mlr.use_ss ==0)
        return;
    use_ss = 1;
    n = mlr.n;
    p = mlr.p;

    if(mlr.R != 0){
        R = gsl_matrix_calloc(p,p);
        gsl_matrix_memcpy(R, mlr.R);
    }
    if(mlr.Z != 0){
        Z = gsl_matrix_calloc(p,1);
        gsl_matrix_memcpy(Z,mlr.Z);
    }
    sigma_option = mlr.sigma_option;
    phi2_vec = mlr.phi2_vec;


}


void MLR::copy_ss(const MLR & mlr, vector<int> & indicator){
    if(mlr.use_ss == 0)
        return;
    use_ss = 1;
    n = mlr.n;

    int tp = mlr.p;
    int ep = 0;
    int count = 0;
    std::map<int,int> imap;


    for(int i=0;i<indicator.size();i++){
        if(indicator[i]==1){
            ep++;
            imap[i] = count++;
        }
    }
    
    p = ep;
    R = gsl_matrix_calloc(ep,ep);
    Z = gsl_matrix_calloc(ep,1);


    for(int i=0;i<tp;i++){
        if(indicator[i] == 0)
            continue;
        gsl_matrix_set(Z, imap[i],0,gsl_matrix_get(mlr.Z,i,0));
        for(int j=0;j<tp;j++){
            if(indicator[j] == 1){

                double val = gsl_matrix_get(mlr.R,i,j);
                gsl_matrix_set(R,imap[i], imap[j],val);
            }
        }
    }
    sigma_option = mlr.sigma_option;
    phi2_vec = mlr.phi2_vec;


}




void MLR::copy_full(const MLR& mlr){

    if(mlr.use_ss == 1)
        return;

    use_ss = 0;	
    n = mlr.n;
    p = mlr.p;
    yty = mlr.yty;

    if(mlr.GtG !=0){
        GtG = gsl_matrix_calloc(p,p);
        gsl_matrix_memcpy(GtG, mlr.GtG);
    }
    if(mlr.Gty!=0){
        Gty = gsl_matrix_calloc(p,1);
        gsl_matrix_memcpy(Gty, mlr.Gty);
    }
    if(mlr.G!=0){
        G = gsl_matrix_calloc(n,p);
        gsl_matrix_memcpy(G,mlr.G);
    }	

    if(mlr.Y!=0){
        Y=gsl_matrix_calloc(n,1);
        gsl_matrix_memcpy(Y,mlr.Y);
    }

    sigma_option = mlr.sigma_option;
    phi2_vec = mlr.phi2_vec;
}

void MLR::copy_full(const MLR & mlr, vector<int> & indicator){
   
    if(mlr.use_ss == 1)
        return;

    use_ss = 0;
    n =  mlr.n;
    yty = mlr.yty;

    int tp = mlr.p;
    int ep = 0;
    int count = 0;
    std::map<int,int> imap;


    for(int i=0;i<indicator.size();i++){
        if(indicator[i]==1){
            ep++;
            imap[i] = count++;
        }
    }
    
    p = ep;

    if(mlr.GtG !=0){
        GtG = gsl_matrix_calloc(p,p);
        for(int i=0;i<tp;i++){
            if(indicator[i] == 0)
                continue;
            for(int j=0;j<tp;j++){
                if(indicator[j] == 1){
                    double val = gsl_matrix_get(mlr.GtG,i,j);
                    gsl_matrix_set(GtG,imap[i], imap[j],val);
                }   
            }   
        }   
    }

    if(mlr.Gty!=0){
        Gty = gsl_matrix_calloc(p,1);
        for(int i=0;i<tp;i++){
            if(indicator[i] == 0)
                continue;
            gsl_matrix_set(Gty, imap[i], 0, gsl_matrix_get(mlr.Gty, i,0));
        }
    }

    if(mlr.G!=0){
        G = gsl_matrix_calloc(n,p);
        for(int i=0;i<tp;i++){
            if(indicator[i] == 0)
                continue;
            for(int j=0;j<n;j++){
                gsl_matrix_set(G, j,imap[i] , gsl_matrix_get(mlr.G, j,i));
            }

        }

    }

    if(mlr.Y!=0){
        Y=gsl_matrix_calloc(n,1);
        gsl_matrix_memcpy(Y,mlr.Y);
    }

    sigma_option = mlr.sigma_option;
    phi2_vec = mlr.phi2_vec;

}



MLR::~MLR(){
    if(GtG!=0)
        gsl_matrix_free(GtG);
    if(Gty!=0)
        gsl_matrix_free(Gty);

    if(Y!=0)
        gsl_matrix_free(Y);
    if(G!=0)
        gsl_matrix_free(G);
    if(Z!=0)
        gsl_matrix_free(Z);
    if(R!=0)
        gsl_matrix_free(R);

}

// ================ Setting parmaeters/options ======================== //


void MLR::set_effect_vec(vector<double> &phi2_vec_){ 

    phi2_vec = phi2_vec_;

}


double MLR::compute_log10_ABF(){
   vector<int> indicator(p,1);
   return(compute_log10_ABF(indicator));
}



double MLR::compute_log10_ABF(vector<int> & indicator){
    if(!use_ss){
        if(sigma_option>=0 && sigma_option <=1)
            return compute_log10_ABF_FD(indicator);
        else
            return compute_log10_BF_FD(indicator);
    }else
        return compute_log10_ABF_SS(indicator);


}





double MLR::compute_log10_BF_FD(vector<int> & indicator){
    vector<double> rstv;
    vector<double> wv;


    vector<int> index_vec;

    for(int i=0;i<indicator.size();i++){
        if(indicator[i]==1){
            index_vec.push_back(i);
        }
    }

    int ep = index_vec.size();

    gsl_matrix *X = gsl_matrix_calloc(n, ep);
    gsl_matrix *XtX = gsl_matrix_calloc(ep,ep);
    gsl_matrix *Xty = gsl_matrix_calloc(ep,1);

    gsl_vector *nv = gsl_vector_alloc(n);

    for(int i=0;i<index_vec.size();i++){
        int index = index_vec[i];
        gsl_matrix_get_col(nv, G, index);
        gsl_matrix_set_col(X, i, nv);
    }


    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0, XtX);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,Y,0, Xty);


    // SV decomp of XtX
    gsl_matrix *V = gsl_matrix_calloc(ep,ep);
    gsl_vector *S = gsl_vector_calloc(ep);
    gsl_vector *work = gsl_vector_calloc(ep);
    gsl_linalg_SV_decomp (XtX, V, S,work);

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

        double b_rss = gsl_matrix_get(tt4,0,0);
        double log10BF = -0.5*log10(fabs(det)) - 0.5*n*log10(1 - b_rss/yty);

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

    gsl_matrix_free(XtX);
    gsl_matrix_free(Xty);

    gsl_vector_free(nv);
    gsl_matrix_free(X);

    double rst =  log10_weighted_sum(rstv,wv);

    return rst;

}


// ABF with full data
double MLR::compute_log10_ABF_FD(vector<int> & indicator){


    // construct the sub-matrices

    vector<double> rstv;
    vector<double> wv;


    vector<int> index_vec;

    for(int i=0;i<indicator.size();i++){
        if(indicator[i]==1){
            index_vec.push_back(i);
        }
    }

    int ep = index_vec.size();

    gsl_matrix *X = gsl_matrix_calloc(n, ep);
    gsl_matrix *XtX = gsl_matrix_calloc(ep,ep);
    gsl_matrix *Xty = gsl_matrix_calloc(ep,1);

    gsl_vector *nv = gsl_vector_alloc(n);

    for(int i=0;i<index_vec.size();i++){
        int index = index_vec[i];
        gsl_matrix_get_col(nv, G, index);
        gsl_matrix_set_col(X, i, nv);
    }


    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,X,0, XtX);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,X,Y,0, Xty);



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

    gsl_matrix_free(X);
    gsl_vector_free(nv);

    double rst =  log10_weighted_sum(rstv,wv);

    return rst;
} 

double MLR::compute_log10_ABF_SS(vector<int> & indicator){

    // construct the sub-matrices

    vector<double> rstv;

    int ep = 0;
    int count = 0;
    std::map<int,int> imap;


    for(int i=0;i<indicator.size();i++){
        if(indicator[i]==1){
            ep++;
            imap[i] = count++;
        }
    }


    gsl_matrix *Rm = gsl_matrix_calloc(ep,ep);
    gsl_matrix *Zm = gsl_matrix_calloc(ep,1);


    for(int i=0;i<p;i++){
        if(indicator[i] == 0)
            continue;
        gsl_matrix_set(Zm, imap[i],0,gsl_matrix_get(Z,i,0));
        for(int j=0;j<p;j++){
            if(indicator[j] == 1){

                double val = gsl_matrix_get(R,i,j);
                gsl_matrix_set(Rm,imap[i], imap[j],val);
            }
        }
    }
    gsl_matrix *V = gsl_matrix_calloc(ep,ep);
    gsl_vector *S = gsl_vector_calloc(ep);
    gsl_vector *work = gsl_vector_calloc(ep);
    gsl_linalg_SV_decomp (Rm, V, S,work);


    for(int i=0;i<phi2_vec.size();i++){
        double kappa = phi2_vec[i];
        gsl_matrix *t1 = gsl_matrix_calloc(ep,ep);
        double log_det = 0;
        for(int i=0;i<ep;i++){
            double v = gsl_vector_get(S,i);
            gsl_matrix_set(t1,i,i,1.0/(v+1/kappa));
            log_det += log(1+kappa*v);
        }

        gsl_matrix *t2 = gsl_matrix_calloc(ep,ep);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);


        gsl_matrix *M_inv = gsl_matrix_calloc(ep,ep);
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,M_inv);

        gsl_matrix *t3 = gsl_matrix_calloc(1,ep);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans,1, Zm,M_inv,0,t3);

        gsl_matrix *t4 = gsl_matrix_calloc(1,1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,t3,Zm,0,t4);


        double log_BF = -0.5*log_det + 0.5*gsl_matrix_get(t4,0,0);
        /*
           double K = n*kappa;
           double zs = gsl_matrix_get(Zm,0,0);
           double log_BF2 = 0.5*log(1.0/(1+K))+0.5*(K/(1+K))*zs*zs;
           printf("%f            %f  %f\n", K, log_BF, log_BF2);
           */
        rstv.push_back(log_BF/log(10));

        gsl_matrix_free(t1);
        gsl_matrix_free(t2);
        gsl_matrix_free(t3);
        gsl_matrix_free(t4);
        gsl_matrix_free(M_inv);
    }

    gsl_matrix_free(V);
    gsl_matrix_free(Rm);
    gsl_matrix_free(Zm);

    gsl_vector_free(S);
    gsl_vector_free(work);

    vector<double> wv(phi2_vec.size(), 1.0/phi2_vec.size());
    double rst =  log10_weighted_sum(rstv,wv);
    return rst;




}




// utility

double MLR::log10_weighted_sum(vector<double> &vec, vector<double> &wts){


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





void MLR::print_matrix(gsl_matrix *M, int a, int b, FILE *outfd){

    if(outfd == 0){
        outfd = stdout;
    }     

    for(int i=0;i<a;i++){
        for(int j=0;j<b;j++){
            //printf("%f  ",gsl_matrix_get(M,i,j));
            fprintf(outfd, "%f",gsl_matrix_get(M,i,j));
            if(j!=b-1)
                fprintf(outfd," ");  
        }
        fprintf(outfd,"\n");
    }
}

void MLR::extract_summary(){
    if(use_ss == 1)
        return;

    // output R matrix
    GtG = gsl_matrix_calloc(p,p);
    Gty = gsl_matrix_calloc(p,1);
    R = gsl_matrix_calloc(p,p);
    
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,G,G,0,GtG);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,G,Y,0,Gty);
    
    for(int i=0;i<p;i++){
        for(int j = 0;j<p;j++){
            double val = gsl_matrix_get(GtG,i,j);
            val = val/sqrt(gsl_matrix_get(GtG,i,i)*gsl_matrix_get(GtG,j,j));
            gsl_matrix_set(R,i,j,val);
        }
    }	  

    double Syy = yty;
    Z = gsl_matrix_calloc(p,1);
    for(int i=0;i<p;i++){
        double Sxx = gsl_matrix_get(GtG, i,i);
        double Sxy = gsl_matrix_get(Gty, i,0);
        double s2 = (Syy - Sxy*Sxy/Sxx)/(n-2);
        double beta = Sxy/Sxx;
        double se = sqrt(s2/Sxx); 
        beta_vec.push_back(beta);
        se_vec.push_back(se);
        gsl_matrix_set(Z,i,0,beta/se);
    }

    gsl_matrix_free(GtG);
    gsl_matrix_free(R);
    gsl_matrix_free(Gty);
    GtG = R = Gty = 0;


}


void MLR::get_single_SNP_stats(){
    if(use_ss==1){
        for( int i=0;i<p;i++){
            beta_vec.push_back(gsl_matrix_get(Z,i,0));
            se_vec.push_back(1.0);
        }
        return;
    }
    //else
    double Syy = yty;
    for(int i=0;i<p;i++){
        double Sxx = 0;
        double Sxy = 0;
        for(int j=0;j<n;j++){
            double gv = gsl_matrix_get(G,j,i);
            Sxx += gv*gv;
            Sxy += gv*gsl_matrix_get(Y,j,0);
        }
        double s2 = (Syy - Sxy*Sxy/Sxx)/(n-2);
        double beta = Sxy/Sxx;
        double se = sqrt(s2/Sxx);
        beta_vec.push_back(beta);
        se_vec.push_back(se);
    }
    return;
}



