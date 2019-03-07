#include "parser.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "MLR.h"
#include <set>
#include <map>
#include <stdio.h>

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


        int max_size; //user option, maximum model size


        // collection of all sizes of size_model
        vector<size_model> szm_vec; // size_model collections


        //candidate set/map of SNPs
        vector<int> cand_set;  // candidate SNP set for higher order models
        map<int, int> cand_map;
        map<int, int> snp2cluster_map;

        // for reporting
        map<string,int> nsnp_map;

        vector<NSNP> nsnp_vec;
        vector<Nmodel> nmodel_vec;


        // output pip for all SNPs default only output those > 0.001
        int output_all;


        // log10 normalizing constant
        double log10_pnorm;


        map<string,double> single_log10_abfv;

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

    public:

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



        void set_prior(char *prior_file);
        void set_prior_exp(double pes);
        void set_prior(double pi1);

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


        // utility
        void parse_nmodel(Nmodel nmod);
        vector<int> get_config(int c);
        double *get_weights(vector<double>& vec);
        double compute_average_r2 (const vector<int> & vec);
        double compute_average_r2(const vector<int> & vec1, const vector<int> & vec2);
        double compute_r2(int i, int j); 
};  

double log10_weighted_sum(vector<double> &vec, vector<double> &wts);

