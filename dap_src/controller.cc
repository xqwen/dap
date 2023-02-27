using namespace std;

#include "controller.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>

void controller::initialize(char *effect_file, char *ld_file, char *grid_file, int sample_size, double syy_, int ld_format){
    
    use_ss = 2;
    N = sample_size;
    syy = syy_;

    if(strlen(grid_file) == 0)
        set_default_grid();
    else
        load_grid(grid_file);
    
    pars.process_summary_data2(effect_file, ld_file,sample_size, syy, ld_format);
    p = pars.ld_matrix->size1;
    set_default_options();
}

void controller::initialize(char *zval_file, char *ld_file, char *grid_file, int sample_size, int ld_format){

    use_ss = 1;
    if(strlen(grid_file) == 0)
        set_default_grid();  
    else
        load_grid(grid_file);
    pars.process_summary_data(zval_file, ld_file, sample_size, ld_format);
    p = pars.ld_matrix->size1;
    N = sample_size;
    // set an approximate smaple size to avoid errors
    if(sample_size <= 0)
        sample_size = 1000;

    set_default_options();
}


void controller::initialize(char *data_file, char *grid_file){

    use_ss = 0;
   
    if(strlen(grid_file) == 0)
        set_default_grid();
    else
        load_grid(grid_file);

    pars.process_data(data_file);

    p = pars.geno_vec[0].size();
    N = pars.geno_vec[0][0].size();   

    // default maximum model size
    set_default_options();
    TWAS_wts = vector<double> (p,0);

}


void controller::set_default_options(){

    max_size = p;
    // default thread 
    nthread = 1;

    use_dap1 = 0;
    output_all = 0;
    size_select_thresh = 0.01;  // set to big positive numbers to enforce the adaptive stopping rule
    //snp_select_thresh = 0.01;


    cluster_pip_thresh = 0.25; // for output cluster purpose
    priority_msize = 100;
    log10_snp_thresh = 2;
    
    // 02/21/2023  set default ld control threshold to 0.5
    ld_control_thresh = 0.5; // by default  ld control with r^2 = 0.5

    //02/26/2023
    csfd = 0;
    coverage_prob = 1.0;

}


// =============== summarize all configuations ======================


void controller::print_dap_config(){
    fprintf(logfd,"\n============ DAP Configuration ============\n\n");
    fprintf(logfd,"INPUT\n\n");
    if(use_ss){
        fprintf(logfd,"\t* summary-level data (z-scores and LD matrix)\n");
        fprintf(logfd,"\t* number of candidate SNPs: %d\n",p);
    }else{
        fprintf(logfd,"\t* individual-level data (sbams format)\n");
        fprintf(logfd,"\t* number of candidate SNPs: %d\n",p);
        fprintf(logfd,"\t* sample size: %d\n",int(pars.geno_vec[0][0].size()));
        fprintf(logfd,"\t* number of covariates: %d\n", int(pars.covar_vec[0].size()));
    }
    fprintf(logfd,"\n");
    
    fprintf(logfd, "PROGRAM OPTIONS\n\n");
    
    if(run_option==0){
        fprintf(logfd,"\t* maximum model size allowed [-msize]: %d ",max_size);
        if(max_size ==p){
            fprintf(logfd,"(no restriction)\n");
        }else{
            fprintf(logfd,"\n");
        }

        fprintf(logfd,"\t* LD control threshold [-ld_control]: %.2f\n", ld_control_thresh);
        fprintf(logfd,"\t* normalizing constant convergence threshold [-converg_thresh]: %.2e (log10 scale)\n", size_select_thresh);
    }

    
   if(run_option==0)
       fprintf(logfd,"\nRUN LOG\n");


}




void controller::load_grid(char *grid_file){
    ifstream gfile(grid_file);
    string line;
    istringstream ins;

    while(getline(gfile,line)){

        ins.clear();
        ins.str(line);
        double omg;
        if(ins>>omg){
            if(use_ss!=1)
                omg2_vec.push_back(pow(omg,2));
            else
                kv_vec.push_back(omg);
        }
    }
    gfile.close();

}

void controller::set_default_grid(){
    if(use_ss!=1){
        omg2_vec.push_back(0.04);
        omg2_vec.push_back(0.16);
        omg2_vec.push_back(0.64);
    }else{
        kv_vec.push_back(4);
        kv_vec.push_back(16);
        kv_vec.push_back(25);
    }


}


void controller::set_outfile(char *outfile, char *logfile){

    if(strlen(outfile)==0)
        outfd = stdout;
    else
        outfd = fopen(outfile,"w");

    
    if(strlen(logfile)==0)
        logfd = stderr;
    else
        logfd = fopen(logfile,"w");
}

void controller::set_cs(char *cs_file, double cs_thresh){
    if(strlen(cs_file)>0)
        csfd = fopen(cs_file, "w");
    coverage_prob = cs_thresh;
}



// exchangeable prior pes/p
void controller::set_prior_exp(double pes){

    double pi1 = pes/double(p);
    set_prior(pi1);
}



// exchangeable prior pi1
void controller::set_prior(double pi1){

    if(pi1 > 1 - 1e-3)
        pi1 = 1-1e-3;

    for(int i=0;i<p;i++){
        pi_vec.push_back(pi1);
    }
    return;
}

void controller::set_prior(char *prior_file){

    if(prior_file==0)
        return;


    ifstream ifile(prior_file);
    string line;
    istringstream ins;

    pi_vec = vector<double>(p,0.0);

    while(getline(ifile,line)){

        ins.clear();
        ins.str(line);
        double value;
        string snp_name;
        ins>>snp_name;

        int index = pars.geno_rmap[snp_name];
        ins>>value;
        pi_vec[index] = value;

    }

    return;


}


void controller::init(){  
    if(use_ss==0){
        vector<vector<vector<double> > > Xgv;

        Xgv = pars.geno_vec;

        if(Xgv.size() == 0)
            exit(0);
        // set SSLR parameters and options
        mlr.set_abf_option(abf_sigma_option); 
        mlr.init(pars.pheno_vec[0],Xgv[0]);

        //sslr.set_prior_option(sslr_prior_option);
        mlr.set_effect_vec(omg2_vec);
        
    }else if(use_ss==1){
        mlr.init_ss(pars.zval_matrix, pars.ld_matrix);
        mlr.set_effect_vec(kv_vec);
    }else if(use_ss==2){
        gsl_matrix *GtG = gsl_matrix_calloc(p,p);
        gsl_matrix_memcpy(GtG,pars.ld_matrix);
        gsl_matrix *Gty = gsl_matrix_calloc(p,1);

        for(int i=0;i<p;i++){
            gsl_matrix_set(Gty, i,0, pars.gty_vec[i]);
        }

        for(int i=0;i<p;i++){
            for(int j=0;j<p;j++){
                double factor = sqrt(pars.sxx_vec[i]*pars.sxx_vec[j]);
                double val = factor*gsl_matrix_get(GtG,i,j);
                gsl_matrix_set(GtG,i,j,val);
            }
        }
        mlr.init(syy, GtG, Gty, N);
        mlr.set_abf_option(abf_sigma_option);
        
        mlr.set_effect_vec(omg2_vec);
    }

    null_config = vector<int>(p,0);  
    //start with empty table

    double sum = 0;
    for(int i=0;i<p;i++){
        sum += (pi_vec[i])/(1-pi_vec[i]);
    }

    prior_ratio = sum/p;

}




double controller::compute_log10_prior(vector<int> &mcfg){

    double lp=0;
    for(int i=0;i<p;i++){
        if(mcfg[i]==0) 	 
            lp += log(1-pi_vec[i]);
        else
            lp += log(pi_vec[i]);
    }

    return lp/log(10);

}


void controller::run(){
    switch(run_option){
        case 0:
            fine_map();
            break;
        case 1:
            scan();
            break;
        case 2:
            extract_ss();
            break;
        case 3:
            extract_ss2();
            break;
        default:
            break;
    }
}

void controller::fine_map(){


    init();

    // initialize snp2cluster_map
    for(int i=0;i<p;i++){
        snp2cluster_map[i] = -1;
    }


    vector<double> log10_pmass_vec;



    double curr_val = compute_log10_prior(null_config);
    log10_pmass_vec.push_back(curr_val);


    fprintf(logfd, "\n\n\tModel_Size \tcandidates \t  log10(NC)\n");
   
    // record best model
    vector<int> bm;
    // single SNP scan, label candidate SNPs for higher order models
    szm_vec.push_back(compute_post_model_single(bm));
    log10_pmass_vec.push_back(szm_vec[0].log10_sum_post);
    vector<double> wv(log10_pmass_vec.size(),1.0);
    double val = log10_weighted_sum(log10_pmass_vec,wv);
    int total_snp = szm_vec[0].mvec.size();
    fprintf(logfd,  "\t %4d \t \t %4d \t \t %9.3f\n",1,total_snp,val);

    double prev_val = val;
    int size = 1;
    int no_bc = 0;
    double increment;

    while(1){

        // if already at allowed max size, quit search
        if(bm.size() == max_size)
            break;
    

        // attempt to add one more variable
        double log10_post = conditional_est(bm);

        // if no more variable to add, then quit search
        if(size == bm.size())
            break;
        else
            size++;

        
        // perform backwards checking, see if can eliminate one variable 
        int elim_index = -1;
        if(size >= 3 && !no_bc){
            elim_index = backward_checking(bm, log10_post);
            if(elim_index !=-1){
                size--;
            }
        }


        fprintf(logfd,"\t %4d \t \t ",size);
        fflush(logfd);

        size_model szm;

        double cps = szm_vec[szm_vec.size()-1].log10_sum_post;
        // append an equal size model to correct for necessary backward checking
        if(elim_index != -1){
            
            map<int, int> black_list;
            for(int i=0;i<szm_vec[elim_index].snp_cluster.size();i++){
                int snp = szm_vec[elim_index].snp_cluster[i];
                black_list[snp] = 100;
            }
            szm = append_post_model(size, black_list); 
            
            szm.snp_cluster = cand_set;
            szm_vec.push_back(szm);
            for(int i=0;i<cand_set.size();i++){
                snp2cluster_map[cand_set[i]] = szm_vec.size();
            }
            log10_pmass_vec.push_back(szm.log10_sum_post);      
            // update log10(NC)
            vector<double> nwv(log10_pmass_vec.size(),1.0);
            val = log10_weighted_sum(log10_pmass_vec,nwv);
            total_snp += cand_set.size();
            fprintf(logfd,  "%4d \t \t %9.3f\t*\n",total_snp,val);
            //  stop continuous backward checking for tiny improvement
            if(val - prev_val < size_select_thresh)
                no_bc = 1;  

        }else{
            //allow backward checking if adding a new variable    
            if(no_bc)
                no_bc = 0;

            // append a size 
            int use_abs_cutoff = 0;
            if(increment > 1 || size == 2)
                use_abs_cutoff = 1;


            szm = compute_post_model(size, use_abs_cutoff);
            szm.snp_cluster = cand_set;
            szm_vec.push_back(szm);
            for(int i=0;i<cand_set.size();i++){
                snp2cluster_map[cand_set[i]] = szm_vec.size();
            }
            log10_pmass_vec.push_back(szm.log10_sum_post);      
            // update log10(NC)
            vector<double> nwv(log10_pmass_vec.size(),1.0);
            val = log10_weighted_sum(log10_pmass_vec,nwv);
            total_snp += cand_set.size();
            fprintf(logfd,  "%4d \t \t %9.3f\n",total_snp,val);
        }



        // stopping rule

        int cs = size;
        // if the prior expectation is high for number of signals, carry on

        // expected contribution of next size partition, assuming model is saturated 
        double project_ratio = (p-cs+1)*prior_ratio/cs;
        if(project_ratio > 1){
            increment = val - prev_val;
            prev_val = val;
            continue;
        }
        // else check stopping point
        double ncps = szm.log10_sum_post;
        double rb = log10(double(p)-cs+1)+log10(prior_ratio) + cps;
        double lb = log10(double(p-2*cs+2)/cs) + log10(prior_ratio) + cps;
        if( ncps <= rb && val - prev_val <= size_select_thresh){
            break;
        }

        increment = val - prev_val;
        prev_val = val;
        //fprintf(stderr, "increment = %.3f\n", increment);

    }
    summarize_approx_posterior();

}




double controller::conditional_est(vector<int>& bm){

    cand_set.clear();

    double max = -9999;

    vector<int>  mcfg = null_config;


    for(int i=0;i<bm.size();i++){
        int index = bm[i];
        mcfg[index] = 1;
    }
    

    vector<double> post_vec;

    int max_id;
    for(int i=0;i<p;i++){

        if(cand_map[i]==1){
            post_vec.push_back(-99999);
            continue;
        }


        vector<int>  cmcfg = mcfg;


        cmcfg[i] = 1;


        double rst =  mlr.compute_log10_ABF(cmcfg) + compute_log10_prior(cmcfg);
        post_vec.push_back(rst);

        if(rst > max){
            max= rst;
            max_id = i;
        }

    }

    double thresh = -99999;


    if(post_vec.size()>size_limit){
        vector<double> post_vec_sort = post_vec;
        sort(post_vec_sort.rbegin(), post_vec_sort.rend());
        thresh = post_vec_sort[size_limit-1];
    }

    int flag = 0;
    if(max > -9999){
        for(int i=0;i<p;i++){
            if( max - post_vec[i]  <= log10_snp_thresh  && post_vec[i]>=thresh){
                if(compute_r2(max_id,i)< ld_control_thresh)
                    continue;
                cand_set.push_back(i);
                cand_map[i]=1;
                flag =1;
            }

        }
    }

    if(flag == 1)
        bm.push_back(max_id);

    return max;
}


// backwards elimination checking
int controller::backward_checking(vector<int>& bm, double log10_post){


    double max = -99999;
    int max_snp = -1;
    int max_id = -1;


    for(int i=0; i<bm.size()-2; i++){
        vector<int> tm = bm;
        tm.erase(tm.begin()+i);

        vector<int> mcfg = null_config;

        for(int j=0;j<tm.size();j++){
            int index = tm[j];
            mcfg[index] = 1;
        }

        double rst =  mlr.compute_log10_ABF(mcfg) + compute_log10_prior(mcfg);
        if(rst > max){
            max= rst;
            max_snp = bm[i];
            max_id=i;
        }

    }

    if(max-log10_post>0){
        bm.erase(bm.begin()+max_id);
        return snp2cluster_map[max_snp];
    }


    return -1;

}

size_model controller::append_post_model(int size, map<int, int> &black_list){

    size_model smod;
    smod.size = size;
    smod.log10_sum_post = 0;

    vector<vector<int> > & model_vec = szm_vec[szm_vec.size()-1].mvec;
    vector<vector<int> > mc_vec;

    map<string, int> saved_model_string_map;

    for(int i=0; i<model_vec.size(); i++){
        for(int j=0;j<cand_set.size();j++){
            vector<int> cm = model_vec[i];

            for(int k=0;k<cm.size();k++){
                if(black_list[cm[k]] == 100){
                    cm.erase(cm.begin()+k);
                    break;
                }
            }
            if(cm.size()==size-1){

                cm.push_back(cand_set[j]);

                stringstream mstream;
                
                for(int k=0;k<cm.size();k++){
                    mstream<< cm[k]<<":";
                }

                string model_string = mstream.str(); 
                if(saved_model_string_map[model_string] != 100){;
                    mc_vec.push_back(cm);
                    saved_model_string_map[model_string] = 100;
                }

            }
        }
    }

    vector<double> post_vec(mc_vec.size());
    int ms = mc_vec.size();

    map<string, double> post_map;
    vector<string> name_vec(ms);

#pragma omp parallel for num_threads(nthread)  
    for(int i=0;i<ms;i++){

        // empty sets to start
        vector<int>  mcfg = null_config;
        vector<int> cm = mc_vec[i];

        string name ="";
        for(int j=0;j<size;j++){
            int index = cm[j];
            name += pars.geno_map[index];
            if(j!=size-1){
                name += "&";
            }

            mcfg[index]=1;
        }



        MLR local_mlr;

        local_mlr.copy(mlr,mcfg);
        double log10_abf  = local_mlr.compute_log10_ABF();
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
#pragma omp critical 
        {
            post_vec[i] = log10_post;
            name_vec[i] = name;
            post_map[name] = log10_post;
        }
    }


    vector<double> rst_post_vec;
    map<string, double> rst_post_map;
    vector<vector<int> > rst_mc_vec;

    vector<double> post_vec_sort = post_vec;
    sort(post_vec_sort.rbegin(), post_vec_sort.rend());

    double cutoff = post_vec_sort[0] - 2;

    for(int i=0;i<post_vec.size();i++){
        if(post_vec[i]>= cutoff){
            rst_post_vec.push_back(post_vec[i]);
            rst_mc_vec.push_back(mc_vec[i]);
            rst_post_map[name_vec[i]] = post_vec[i];
        }
    }

    vector<double> wv(rst_post_vec.size(),1.0);
    smod.log10_sum_post = mlr.log10_weighted_sum(rst_post_vec,wv);
    smod.post_map = rst_post_map;
    smod.mvec = rst_mc_vec;

    return smod;




}

// single SNP model, required to run 
size_model controller::compute_post_model_single(vector<int>& bm){



    size_model smod;
    smod.size = 1;
    smod.log10_sum_post = 0;

    vector<double> post_vec;
    vector<double> abf_vec;

    double max_log10_abf = -9999;
    double max_log10_post = -9999;
    int max_id = -1;
    mlr.get_single_SNP_stats();
    for(int index=0;index<p;index++){

        string name = pars.geno_map[index];
        vector<int>  mcfg = null_config;
        mcfg[index]=1;
        double log10_abf = mlr.compute_log10_ABF(mcfg); 
        single_log10_abfv[name] = log10_abf;
        single_bhat[name] = mlr.beta_vec[index];
        single_se[name] = mlr.se_vec[index];
        abf_vec.push_back(log10_abf);
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
        post_vec.push_back(log10_post);

        if(log10_post>max_log10_post){
            max_log10_post = log10_post;
            max_id = index;
        }

        smod.post_map[name] = log10_post;
    }


    vector<double> wv(post_vec.size(),1.0);
    smod.log10_sum_post = mlr.log10_weighted_sum(post_vec,wv);

    double thresh = -99999;

    if(post_vec.size()>size_limit){
        vector<double> post_vec_sort = post_vec;
        sort(post_vec_sort.rbegin(), post_vec_sort.rend());
        thresh= post_vec_sort[size_limit-1];
    }


    for(int i=0;i<p;i++){
        if( (max_log10_post - post_vec[i]  <= log10_snp_thresh && post_vec[i]>= thresh) || use_dap1 ) {
            if(compute_r2(max_id,i)< ld_control_thresh)
                continue;
            vector<int> mv;
            mv.push_back(i);
            smod.mvec.push_back(mv);
            cand_set.push_back(i);
            cand_map[i] = 1;
            snp2cluster_map[i]=0;
        }
    }
    

    smod.snp_cluster = cand_set;
    bm.push_back(max_id);
    return smod;

}


size_model controller::compute_post_model(int size, int use_abs_cutoff){


    size_model smod;
    smod.size = size;
    smod.log10_sum_post = 0; 


    vector<vector<int> > & model_vec = szm_vec[szm_vec.size()-1].mvec;
    vector<vector<int> > mc_vec;

    for(int i=0; i<model_vec.size(); i++){
        for(int j=0;j<cand_set.size();j++){
            vector<int> cm = model_vec[i];
            cm.push_back(cand_set[j]);
            mc_vec.push_back(cm);

        }
    }



    vector<double> post_vec(mc_vec.size());
    int ms = mc_vec.size();



    map<string, double> post_map;
    vector<string> name_vec(ms);

#pragma omp parallel for num_threads(nthread)  
    for(int i=0;i<ms;i++){

        // empty sets to start
        vector<int>  mcfg = null_config;
        vector<int> cm = mc_vec[i];

        string name ="";
        for(int j=0;j<size;j++){
            int index = cm[j];
            name += pars.geno_map[index];
            if(j!=size-1){
                name += "&";
            }

            mcfg[index]=1;
        }



        MLR local_mlr;

        local_mlr.copy(mlr,mcfg);
        double log10_abf  = local_mlr.compute_log10_ABF();
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
#pragma omp critical 
        {
            post_vec[i] = log10_post;
            name_vec[i] = name;
            post_map[name] = log10_post;
        }
    }


    vector<double> rst_post_vec;
    map<string, double> rst_post_map;
    vector<vector<int> > rst_mc_vec;

    vector<double> post_vec_sort = post_vec;
    sort(post_vec_sort.rbegin(), post_vec_sort.rend());

    double cutoff = post_vec_sort[0] - 2;
    if(!use_abs_cutoff){
        if(post_vec.size()>priority_msize){
            cutoff = post_vec_sort[priority_msize-1];
        }else{
            cutoff = post_vec_sort[post_vec_sort.size()-1];
        }
    }

    for(int i=0;i<post_vec.size();i++){
        if(post_vec[i]>= cutoff){
            rst_post_vec.push_back(post_vec[i]);
            rst_mc_vec.push_back(mc_vec[i]);
            rst_post_map[name_vec[i]] = post_vec[i];
        }
    }

    vector<double> wv(rst_post_vec.size(),1.0);
    smod.log10_sum_post = mlr.log10_weighted_sum(rst_post_vec,wv);
    smod.post_map = rst_post_map;
    smod.mvec = rst_mc_vec;
    return smod;

}



///////////////////////////// Summarizing and Reporting Results ///////////////////////////////////


void controller::summarize_approx_posterior(){

    vector<double>  log10_pmass_vec;
    log10_pmass_vec.push_back(compute_log10_prior(null_config));


    for(int i=0;i<szm_vec.size();i++){
        log10_pmass_vec.push_back(szm_vec[i].log10_sum_post);
    }

    vector<double> wv1(log10_pmass_vec.size(),1.0);
    double log10_pip_NC = mlr.log10_weighted_sum(log10_pmass_vec,wv1);

    double val = log10_pmass_vec[log10_pmass_vec.size()-1];

    double sum = 0;
    double ratio = 1;
    int max_model_size = szm_vec[szm_vec.size()-1].size;
    if(max_size == p){
        for(int k=max_model_size+1;k<=p;k++){
            ratio *= (p-k+1)*prior_ratio/k;
            sum += ratio;
        }

        log10_pmass_vec.push_back(val+log10(sum));
    }

    vector<double> wv(log10_pmass_vec.size(),1.0);
    log10_pnorm  = mlr.log10_weighted_sum(log10_pmass_vec,wv);


    Nmodel nm;
    nm.id = "NULL";
    nm.prob = pow(10, log10_pmass_vec[0]-log10_pnorm);
    double null_prob = nm.prob;
    nm.size = 0;
    nm.post_score = log10_pmass_vec[0];
    nmodel_vec.push_back(nm);

    for(int i=0;i<szm_vec.size();i++){
        
        map<string, double>::iterator iter;
        for(iter = szm_vec[i].post_map.begin(); iter != szm_vec[i].post_map.end(); iter++){
            Nmodel nm;
            nm.id = iter->first;
            nm.post_score = iter->second;
            nm.prob = pow(10, iter->second-log10_pnorm);
            nm.size = szm_vec[i].size;
            nmodel_vec.push_back(nm);
            parse_nmodel(nm);
        }
    }

    std::sort(nmodel_vec.begin(),nmodel_vec.end(),sort_nmodel_dec);
    double cump = 0;

    double msize_mean = 0;
    double msize_var = 0;

     fprintf(outfd, "model_rank   posterior_prob   size   log10_posterior_score   model_composition\n");

    for(int i=0;i<nmodel_vec.size();i++){

        string name = nmodel_vec[i].id;

        size_t start_pos = 0;
        while((start_pos = name.find("&", start_pos)) != std::string::npos) {
            name.replace(start_pos, 1 , "] [");
            start_pos += 3; 
        }

        if( nmodel_vec[i].prob >= 1e-5 || name == "NULL"){
            fprintf(outfd, "%5d   %7.4e    %d    %7.3f   [%s]\n",i+1, nmodel_vec[i].prob, nmodel_vec[i].size, nmodel_vec[i].post_score, name.c_str());
        }
        cump += nmodel_vec[i].prob;
        msize_mean += nmodel_vec[i].prob*nmodel_vec[i].size;
        msize_var  += nmodel_vec[i].prob*pow(nmodel_vec[i].size,2);

    }



    msize_var -= pow(msize_mean,2.0);
    if(msize_var < 0){
        msize_var = 0;
    }


    fprintf(outfd, "\nPosterior expected model size: %.3f (sd = %.3f)\n", msize_mean, sqrt(msize_var));
    fprintf(outfd, "LogNC = %.5f ( Log10NC = %.3f )", log10_pnorm/log10(exp(1)), log10_pnorm);


    vector<NSNP> nsnp_vec_sort = nsnp_vec;
    std::sort(nsnp_vec_sort.begin(),nsnp_vec_sort.end(),sort_nsnp_dec_by_ip);
    
    

    map<string, int> snp2index;
    for(int i=0;i<nsnp_vec_sort.size();i++){
        //if(nsnp_vec_sort[i].incl_prob<1e-3&&!output_all)
        //    break;
        nsnp_vec_sort[i].cluster = -1;
        snp2index[nsnp_vec_sort[i].name] = i;
    }
    // estimate min_pip from BIC approximation
    int NN = N;
    if(N<0)
        NN=1000;
    double min_pip = (1-null_prob)*prior_ratio/sqrt(NN);
    fprintf(outfd, "\nMinimum PIP is estimated at %7.3e (N = %d)", min_pip, NN);
    if(N<0)
        fprintf(outfd,"   *** assumed sample size");
    fprintf(outfd, "\n");

    vector<double> cluster_pip;
    vector<double> cluster_r2;
    vector<int> cluster_count;
    vector<int> cluster_id;
    vector<vector<int> > grp_vec;
    map<int,int> grpr2_map;
    int cluster_index = 1;

    // organize SNPs by cluster id
    map<int, vector<NSNP> > cluster_map;

    for(int i=0;i<szm_vec.size();i++){

        double cluster_prob = 0;
        vector<int> member_vec;

        for(int j=0;j<szm_vec[i].snp_cluster.size();j++){
            int snp = szm_vec[i].snp_cluster[j];
            string sname = pars.geno_map[szm_vec[i].snp_cluster[j]];
            if(snp2index.find(sname) == snp2index.end())
                continue;
            int index = snp2index[sname];
            double prob = nsnp_vec_sort[index].incl_prob;  
            // we don't consider a SNP as a noteworthy cluster memeber if its pip < min_pip
            
            if(prob > min_pip || use_dap1 ){
                member_vec.push_back(snp);
                cluster_prob += prob;
                nsnp_vec_sort[index].cluster = cluster_index;
            }
            
        }
        if(member_vec.size()>0){
            
            cluster_r2.push_back(compute_average_r2(member_vec));
            cluster_count.push_back(member_vec.size());
            cluster_pip.push_back(cluster_prob);
            cluster_id.push_back(cluster_index);

            if(cluster_prob >= cluster_pip_thresh){
                grp_vec.push_back(member_vec);
                grpr2_map[cluster_count.size()-1] = grp_vec.size()-1;
            }
            cluster_index++;
        }

    }       

    map<string, double> grp_r2;
    if(grp_vec.size()>1){

        for(int i=0;i<grp_vec.size();i++){
            for(int j=i+1;j<grp_vec.size();j++){
                string id(i +":"+j);
                double r2 = compute_average_r2(grp_vec[i], grp_vec[j]);
                grp_r2[id] = r2;
            }
        }
    }

    fprintf(outfd,"\nPosterior inclusion probability\n\n");
    fprintf(outfd,"variant_rank\tvariant_id\tpip\tlog10_BF(multi_SNP)\tsignal_cluster_id\tTWAS_weight\n");

    for(int i=0;i<nsnp_vec_sort.size();i++){
        if(nsnp_vec_sort[i].incl_prob < min_pip && !use_dap1)
            nsnp_vec_sort[i].incl_prob = min_pip;

        string snp_id = nsnp_vec_sort[i].name;
        int cluster_id = nsnp_vec_sort[i].cluster;
        double pip =  nsnp_vec_sort[i].incl_prob;
        if(pip>1-1e-8)
            pip = 1 - 1e-8;
        int index = pars.geno_rmap[nsnp_vec_sort[i].name];
        double weight =TWAS_wts[index];
        if(fabs(weight) < 1e-6)
            weight = 0;
        
        double log10_BF =  log10(pip) - log10(1-pip) - (log10(pi_vec[index]) -log10(1-pi_vec[index]) );
        //double log10_BF = single_log10_abfv[nsnp_vec_sort[i].name]; // alternative to output single SNP BF

        if(cluster_id != -1){
            if(cluster_map.find(cluster_id) == cluster_map.end() ){
                cluster_map[cluster_id] = vector<NSNP>();
            }
            NSNP nsnp;
            nsnp.name = snp_id;
            nsnp.incl_prob = pip;
            cluster_map[cluster_id].push_back(nsnp);
        } 

        // since version 2, the default is output_all, the next two conditions are ineffective
        if(nsnp_vec_sort[i].incl_prob<1e-3&&!output_all)
            break;
        if(nsnp_vec_sort[i].cluster==-1&&!output_all)
            continue;
        /////////////////////////////////////////////////////////////////////////////////////    

        //fprintf(outfd,"((%d))\t%15s\t%8.5e\t%7.3f\t\t%d\t\t%7.3e\n",i+1, nsnp_vec_sort[i].name.c_str(), pip, , nsnp_vec_sort[i].cluster, weight);
        fprintf(outfd,"((%d))\t%15s\t\t%8.3e\t%7.3f\t\t%d\t\t%11.3e\n",i+1, nsnp_vec_sort[i].name.c_str(), pip, log10_BF, cluster_id, weight);
    }


    if(grp_vec.size()>0){
        fprintf(outfd, "\nIndependent association signal clusters\n\n");
        fprintf(outfd,"\t cluster \t member_snp\t cluster_pip\t  average_r2\n");


        for(int i=0;i<cluster_count.size();i++){
            if(cluster_pip[i]<cluster_pip_thresh)
                continue;
            fprintf(outfd,"\t   {%d}\t\t  %3d\t\t %7.3e \t  %7.3f    \t\t",i+1,cluster_count[i],cluster_pip[i],cluster_r2[i]);
            int counti = grpr2_map[i];

            for(int k=0;k<cluster_count.size();k++){
                if(grpr2_map.find(k) == grpr2_map.end())
                    continue;
                int countk = grpr2_map[k];
                string id; 
                if(countk == counti){
                    fprintf(outfd,"%7.3f\t",cluster_r2[counti]);
                    continue;
                }
                if(countk<counti)
                    id = countk+ ":"+counti;
                if(countk>counti)
                    id = counti+":"+countk;
                fprintf(outfd,"%7.3f\t",grp_r2[id]);

            }
            fprintf(outfd,"\n");
            counti++;
        }
    }

    if(csfd!=0){
        
        //fprintf(csfd, "\n\n\ncluster_id\t\tSNP\t\tPIP\n");
        for (int i=0;i<cluster_count.size();i++){
            if(cluster_pip[i] < coverage_prob )
              continue;
            
            int cid = i+1;
            double sum = 0;

            std::sort(cluster_map[cid].begin(),cluster_map[cid].end(),sort_nsnp_dec_by_ip);
            for(int k=0;k<cluster_map[cid].size();k++){
                string snp_id = cluster_map[cid][k].name;
                double pip = cluster_map[cid][k].incl_prob;
                fprintf(csfd,"%2d\t\t%15s\t\t%8.3e\n",cid,snp_id.c_str(),pip);
                sum += pip;
                if(sum>coverage_prob)
                    break;
            }
        }

    
        
        //double log10_BF = single_log10_abfv[nsnp_vec_sort[i].name];
    }





}



void controller::parse_nmodel(Nmodel nmod){

    istringstream iss(nmod.id);
    string token;
    vector<int> index_vec;
    while (getline(iss, token, '&')){
        string snp_id = token;
        int index;
        if(nsnp_map.find(snp_id)!=nsnp_map.end()){
            index = nsnp_map[snp_id];
        }else{
            NSNP ns;
            ns.name = snp_id;
            ns.incl_prob = 0;
            nsnp_vec.push_back(ns);
            index = nsnp_vec.size()-1;
            nsnp_map[snp_id] = index;
        }
        index_vec.push_back(pars.geno_rmap[snp_id]);
        nsnp_vec[index].incl_prob += nmod.prob;
    }
    vector<double> coef_vec = mlr.fit_LS(index_vec);
    for(int i=0;i<index_vec.size();i++){
        TWAS_wts[index_vec[i]] += nmod.prob*coef_vec[i];
    }


}








double *controller::get_weights(vector<double>& vec){

    double max = vec[0];
    for(int i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(int i=0;i<vec.size();i++){
        sum += pow(10, (vec[i]-max));
    }
    double *pp = new double[vec.size()];
    for(int i=0;i<vec.size();i++){
        pp[i] = pow(10, (vec[i]-max))/sum;
    }
    return pp;

}


// utility for single SNP analysis -- input individual-level data output: log10 BF
void controller::scan(){

    init();
    mlr.get_single_SNP_stats();
    for(int i=0;i<p;i++){

        vector<int>  mcfg = null_config;
        mcfg[i] = 1;
        double bhat = mlr.beta_vec[i];
        double se = mlr.se_vec[i];
        double zval = bhat/se;
        fprintf(outfd,"%25s  %9.3f        %12.3e %12.3e   %12.3e\n" , pars.geno_map[i].c_str(),  mlr.compute_log10_ABF(mcfg),bhat, se, zval);

    }  

}


void controller::extract_ss(){

    if(use_ss){
        fprintf(stderr, "Error: can't extract summary info without individual-level data\n");
        exit(1);
    }

    init();
    string z_file = "zval.dat";
    string R_file = "LD.dat";
    if(gene !=""){
        z_file = gene+"."+z_file;
        R_file = gene+"."+R_file;
    }
    mlr.extract_summary();
    gsl_matrix* R = mlr.get_R();
    gsl_matrix* Z = mlr.get_Z();
    FILE *fd = fopen(R_file.c_str(),"w");
    mlr.print_matrix(R, p,p, fd);
    fclose(fd);

    fd = fopen(z_file.c_str(),"w");
    for(int i=0;i<p;i++){
        string name = pars.geno_map[i];
        //fprintf(fd,"%15s\t%7.3e\n",name.c_str(), gsl_matrix_get(Z,i,0));
        fprintf(fd,"%s %f\n",name.c_str(), gsl_matrix_get(Z,i,0));

    }
    fclose(fd);

}

void controller::extract_ss2(){

    if(use_ss){
        fprintf(stderr, "Error: can't extract summary info without individual-level data\n");
        exit(1);
    }

    init();
    string est_file = "est.dat";
    string R_file = "LD.dat";

    if(gene !=""){
        est_file = gene+"."+est_file;
        R_file = gene+"."+R_file;
    }
    mlr.extract_summary();
    gsl_matrix * R = mlr.get_R();
    FILE *fd = fopen(R_file.c_str(),"w");
    mlr.print_matrix(R, p,p, fd);
    fclose(fd);
    fd = fopen(est_file.c_str(),"w");
    for(int i=0;i<p;i++){
        string name = pars.geno_map[i];
        fprintf(fd,"%15s\t%7.3e\t%7.3e\n",name.c_str(), mlr.beta_vec[i], mlr.se_vec[i]);
    }
    fclose(fd);
    mlr.get_summary();
}


double controller::compute_average_r2(const vector<int> & vec){
    if(vec.size()==1)
        return 1.0;
    double sum_r2 = 0;
    int count = 0;
    for(int i=0;i<vec.size();i++){
        for(int j=i+1;j<vec.size();j++){
            sum_r2 += compute_r2(vec[i],vec[j]);
            count++;
        }
    }
    return sum_r2/count;
}


double controller::compute_average_r2(const vector<int> & vec1, const vector<int> & vec2){

    double sum_r2 = 0;
    int count = 0;
    for(int i=0;i<vec1.size();i++){
        for(int j=0;j<vec2.size();j++){
            sum_r2 += compute_r2(vec1[i],vec2[j]);
            count++;
        }
    }
    return sum_r2/count;
}




double controller::compute_r2(int i, int j){
    double r2;
    if(use_ss==0){
        double *gi = &pars.geno_vec[0][i][0];
        double *gj = &pars.geno_vec[0][j][0];
        r2 = pow(gsl_stats_correlation(gi,1,gj, 1, pars.geno_vec[0][i].size()),2);
    }else{
        r2 = pow(gsl_matrix_get(pars.ld_matrix,i,j),2);
    }
    return r2;
}






double log10_weighted_sum(vector<double> &vec, vector<double> &wts){


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









bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs){
    return lhs.incl_prob > rhs.incl_prob;
}




bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs){
    return lhs.prob > rhs.prob;
}


