#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <map>


class MLR {

    private:


        int p; // genotype size
        int n; // sample size

        // full summary data input
        gsl_matrix *GtG;
        gsl_matrix *Gty;
        double yty;

        // full data input
        gsl_matrix *Y;
        gsl_matrix *G;


        // partial summary data input
        gsl_matrix *Z;
        gsl_matrix *R;



        vector<double> phi2_vec; //effect size grid full data
        vector<double> kv_vec; //effect size grid summary-level data

        int use_ss;

    private:
        // options
        double sigma_option;  // 0 to 1, mixture fraction of mle of Sigma under the alternative model, default 0



    public:

        // interface
        // empty constructor, assign default options
        MLR(){
            sigma_option = 0.0; 
            GtG= Gty = 0;
            G = Y = 0;
            Z = R = 0;
            use_ss = 0;
        }

        ~MLR();

        // Note: assuming covariates are regressed beforehand, typically in parser.
        // init
        void init(double yty_, const gsl_matrix *GtG_, const gsl_matrix *Gty_, int n_);
        void init(vector<double> &Y_in, vector<vector<double> >  &G_in);  
        void init_ss(const gsl_matrix *Z_, const gsl_matrix *R_);  
        void copy(const MLR &mlr);

        void copy_full(const MLR &mlr);
        void copy_ss(const MLR &mlr);

        // options
        void set_effect_vec(vector<double> &phi2);
        void set_abf_option(double option){
            sigma_option = option;
        }


        double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

        double compute_log10_ABF(vector<int> &indicator);
        double compute_log10_ABF_FD(vector<int> &indicator);
        double compute_log10_ABF_SS(vector<int> &indicator);

        void extract_summary();
        gsl_matrix * get_R(){
            return R;
        }
        gsl_matrix * get_Z(){
            return Z;
        }

        void print_matrix(gsl_matrix *M, int a, int b, FILE *outfd = 0);
};


