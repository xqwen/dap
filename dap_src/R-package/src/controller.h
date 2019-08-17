#include "parser.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "MLR.h"
#include <set>
#include <map>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

class NSNP {

    public:

        string name;
        double incl_prob;
        int cluster; //signal cluster membership
};

bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs);



class Nmodel {

    public: 
        string id;
        double prob;
        double post_score;
        int size;
};
bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs);



// set of models with same size
class size_model {

    public:
        int size;
        double log10_sum_post;
        map<string, double> post_map;
        vector<vector<int> > mvec; // models to be expanded for the next round
        vector<int> snp_cluster;
};



class controller {

    private:

        FILE *outfd;
        FILE *logfd;

        parser pars;
        string gene;

        vector<double> omg2_vec;
        vector<double> kv_vec;
        MLR mlr;

        int p;

        int N;
        double syy;

        vector<int>  null_config;

        vector<double>  pi_vec; //prior for configs
        double prior_ratio; // for estimate residual from the approximaiton
        double ens; // pes or ens - expected number of signals


        int max_size; //user option, maximum model size


        // collection of all sizes of size_model
        vector<size_model> szm_vec; // size_model collections


        //candidate set/map of SNPs
        vector<int> cand_set;  // candidate SNP set for higher order models
        map<int, int> cand_map;
        map<int, int> snp2cluster_map;

        // for reporting
        map<string,int> nsnp_map;




        // output pip for all SNPs default only output those > 0.001
        int output_all;


        // log10 normalizing constant
        double log10_pnorm;




        // threshold 
        double snp_select_thresh;   // conditional inclusion prob.
        double size_select_thresh;  // the decay of the posteriors of a given size

        // threshold for greedy algorithm
        int    size_limit;
        double cluster_pip_thresh; 
        int    priority_msize;
        double log10_bf_thresh;
        double ld_control_thresh;


        // ABF options
        double abf_sigma_option;

        // openmp thread
        int nthread;
        int use_ss;

        // run in R
        int run_in_r;

    public:
        // for reporting option 0
        vector<Nmodel> nmodel_vec;
        vector<NSNP> nsnp_vec;
        double msize_mean;
        double msize_var;
        double min_pip;

        vector<double> cluster_pip;
        vector<double> cluster_r2;
        vector<int> cluster_count;
        vector<int> cluster_id;
        vector<vector<int> > grp_vec;
        map<int,int> grpr2_map;
        map<string, double> grp_r2;
        map<string,double> single_log10_abfv;

        // for reporting option 1 --scan
        vector<string> single_snp_name;
        vector<double> single_snp_bhat;
        vector<double> single_snp_se;
        vector<double> single_snp_zval;
        vector<double> single_snp_log10_ABF;
        
        // for reporting ld matrix
        vector<double> ld_matrix_flat; // flat vector of size pxp


        // public interface
        int run_option; 
        void initialize(char *effect_file, char *ld_file, char *grid_file, int sample_size, double syy_, int ld_format);
        void initialize(char *data_file, char *grid_file);
        void initialize(char *zval_file, char *ld_file, char *grid_file,int sample_size, int ld_format);

        void run();

        void fine_map();
        void scan();
        void extract_ss();
        void extract_ss2();
        void extract_ss2_in_r();
        void print_dap_config();
        // system setting

        // options

        void set_abf_option(double sigma_option){
            abf_sigma_option= sigma_option;
        }

        void set_size_limit(int size_limit_thresh){
            if(size_limit_thresh<0)
                size_limit = p;
            else
                size_limit = size_limit_thresh;
        }

        void set_ld_control(double thresh){
            ld_control_thresh = thresh;
        }


        void set_ens(double pes){
            ens = pes;
        }
        void set_prior(char *prior_file);
        void set_prior_default();
        void set_prior(double pi1);
        // for R use
        void set_prior(vector<string>& snp_names, vector<double> prior_values);

        void set_max_size(int msize){
            max_size = msize;
        }

        void set_snp_select_thresh (double thresh){
            snp_select_thresh = thresh;
        }

        void set_size_select_thresh (double thresh){
            size_select_thresh = thresh;
        }


        void set_outfile(char *outfile, char *logfile);

        void set_gene(string gname){
            gene=gname;
        }

        void set_output_all(){
            output_all = 1;
            cluster_pip_thresh = -1;

        }


        void set_thread(int thread){
            nthread = thread;
            #ifdef _OPENMP
                omp_set_num_threads(nthread);
            #else
                nthread=1;
            #endif
        }


        // for R use
        void initialize(vector<double>& pheno, string phenoname, vector<vector<double> >& geno, vector<string>& geno_name, char *grid_file, bool regress);
        void initialize(vector<string>& snp, vector<double>& beta, vector<double>& se, vector<vector<double> >& ld, int sample_size, double syy_, string phenoname);
        void initialize(vector<string>& snp, vector<double>& zval, vector<vector<double> >& ld, int sample_size, string phenoname);

        void set_for_r(int quiet){
            //1:not quiet 2:quiet
            run_in_r = 1 + quiet;
        }

        double get_msize_mean(){
            return msize_mean;
        }

        double get_msize_sd(){
            return sqrt(msize_var);
        }

        double get_log10_pnorm(){
            return log10_pnorm;
        }

        double get_min_pip(){
            return min_pip;
        }

        int get_N(){
            return N;
        }

        int get_p(){
            return p;
        }
        
        double get_syy(){
            return syy;
        }

        int get_output_all(){
            return output_all;
        }

        int get_max_size(){
            return max_size;
        }

        double get_ld_control_thresh(){
            return ld_control_thresh;
        }

        double get_size_select_thresh(){
            return size_select_thresh;
        }

        int get_nthread(){
            return nthread;
        }

        string get_pheno_name(){
            return pars.pheno_name;
        }




    private:

        // initialization
        void set_default_options();
        void load_grid(char *grid_file);
        void set_default_grid();

        void init();

        // computing engine
        double compute_log10_prior(vector<int>& mcfg);
        size_model compute_post_model_single(vector<int> &bm);
        size_model compute_post_model(int size, int use_abs_cutoff);
        double conditional_est(vector<int> &control_set);

        int  backward_checking(vector<int>& bm, double log10_post);
        size_model append_post_model(int size, map<int, int> &black_list);

        // reporting
        void summarize_approx_posterior();
        void summary_output();


        // utility
        void parse_nmodel(Nmodel nmod);
        vector<int> get_config(int c);
        double *get_weights(vector<double>& vec);
        double compute_average_r2 (const vector<int> & vec);
        double compute_average_r2(const vector<int> & vec1, const vector<int> & vec2);
        double compute_r2(int i, int j); 

        // for R use
        void init_scan_for_R();
        void test_scan_for_R();
};  

double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

