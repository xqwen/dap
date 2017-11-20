using namespace std;

#include "controller.h"
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>


void controller::initialize(char *zval_file, char *ld_file, char *grid_file){

    use_ss = 1;
    
    if(strlen(grid_file) == 0)
        set_default_grid();  
    else
        load_grid(grid_file);
    pars.process_summary_data(zval_file, ld_file);
    p = pars.ld_matrix->size1;
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

    // default maximum model size
    set_default_options();

}


void controller::set_default_options(){

    max_size = p;
    // default thread 
    nthread = 1;

    output_all = 0;
    size_select_thresh = 0.01;  // set to big positive numbers to enforce the adaptive stopping rule
    snp_select_thresh = 0.01;


    cluster_pip_thresh = 0.25; // for output cluster purpose
    priority_msize = 100;
    log10_bf_thresh = 2;

    ld_control_thresh = 0.25; // by default  ld control with r^2 = 0.25
}

void controller::print_dap_config(){
    
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
            omg2_vec.push_back(pow(omg,2));
        }
    }
    gfile.close();

}

void controller::set_default_grid(){
    if(!use_ss){
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
    if(!use_ss){
        vector<vector<vector<double> > > Xgv;

        Xgv = pars.geno_vec;

        if(Xgv.size() == 0)
            exit(0);
        // set SSLR parameters and options
        mlr.set_abf_option(abf_sigma_option); 
        mlr.init(pars.pheno_vec[0],Xgv[0]);

        //sslr.set_prior_option(sslr_prior_option);
        mlr.set_effect_vec(omg2_vec);
        
    }else{
        mlr.init_ss(pars.zval_matrix, pars.ld_matrix);
        mlr.set_effect_vec(kv_vec);
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


    init();

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

    double increment;

    while(1){
        double log10_abf = conditional_est(bm);
        /*
           printf("%7.3f: [", log10_abf);
           for(int i=0;i<bm.size();i++){
           printf("%d::",bm[i]);
           }
           printf("]\n");
           */
        if(size == bm.size())
            break;
        size++;

        if(size > max_size)
            break;

        fprintf(logfd,"\t %4d \t \t ",size);
        fflush(logfd);

        double cps = szm_vec[szm_vec.size()-1].log10_sum_post;
        int use_abs_cutoff = 0;
        if(increment > 1 || size == 2)
            use_abs_cutoff = 1;

        size_model szm = compute_post_model(size, use_abs_cutoff);
        szm.snp_cluster = cand_set;
        szm_vec.push_back(szm);

        log10_pmass_vec.push_back(szm.log10_sum_post);
        // update log10(NC)
        vector<double> nwv(log10_pmass_vec.size(),1.0);
        val = log10_weighted_sum(log10_pmass_vec,nwv);
        total_snp += cand_set.size();
        fprintf(logfd,  "%4d \t \t %9.3f\n",total_snp,val);




        // stopping rule

        int cs = size;
        // if the prior expectation is high for number of signals, carry on

        // expected contribution of next size partition, assuming model is saturated 
        double project_ratio = (p-cs+1)*prior_ratio/cs;
        if(project_ratio > 1)
            continue;

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
    //vector<double> abf_vec;

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
            if( max - post_vec[i]  <= log10_bf_thresh  && post_vec[i]>=thresh){
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
    for(int index=0;index<p;index++){

        string name = pars.geno_map[index];
        vector<int>  mcfg = null_config;
        mcfg[index]=1;
        double log10_abf = mlr.compute_log10_ABF(mcfg); 
        single_log10_abfv[name] = log10_abf;
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
        if( max_log10_post - post_vec[i]  <= log10_bf_thresh && post_vec[i]>= thresh){
            if(compute_r2(max_id,i)< ld_control_thresh)
                continue;
            vector<int> mv;
            mv.push_back(i);
            smod.mvec.push_back(mv);
            cand_set.push_back(i);
            cand_map[i] = 1;

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

        local_mlr.copy(mlr);
        double log10_abf  = local_mlr.compute_log10_ABF(mcfg);
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
#pragma omp critical 
        {
            post_vec[i] = log10_post;
            name_vec[i] = name;
            post_map[name] = log10_post;
            //printf("\n eval %d  %s  %7.3f %7.3f %7.3f\n", i, name.c_str(), log10_abf, compute_log10_prior(mcfg_map),log10_post);   
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

    if(max_size == p){
        for(int k=log10_pmass_vec.size()+1;k<=p;k++){
            ratio *= (p-k+1)*prior_ratio/k;
            sum += ratio;
        }

        log10_pmass_vec.push_back(val+log10(sum));
    }

    vector<double> wv(log10_pmass_vec.size(),1.0);
    log10_pnorm  = mlr.log10_weighted_sum(log10_pmass_vec,wv);

    /*
       printf("log10_pnorm = %7.3f\n", log10_pnorm);
       for(int i=0;i<wv.size();i++){
       printf("%d:  %7.3f\n", i, log10_pmass_vec[i]);
       }
       */


    Nmodel nm;
    nm.id = "NULL";
    nm.prob = pow(10, log10_pmass_vec[0]-log10_pnorm);
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
            nm.size = i+1;
            nmodel_vec.push_back(nm);
            parse_nmodel(nm);
        }
    }

    std::sort(nmodel_vec.begin(),nmodel_vec.end(),sort_nmodel_dec);
    double cump = 0;

    double msize_mean = 0;
    double msize_var = 0;


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
        if(nsnp_vec_sort[i].incl_prob<1e-3&&!output_all)
            break;
        nsnp_vec_sort[i].cluster = -1;
        snp2index[nsnp_vec_sort[i].name] = i;
    }



    vector<double> cluster_pip;
    vector<double> cluster_r2;
    vector<int> cluster_count;
    vector<vector<int> > grp_vec;
    map<int,int> grpr2_map;
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
            member_vec.push_back(snp);
            nsnp_vec_sort[index].cluster = i+1;
            //printf ("%2d   %15s  %8.5e\n",i+1,sname.c_str(), prob);
            cluster_prob += prob;      
        }
        cluster_r2.push_back(compute_average_r2(member_vec));
        cluster_count.push_back(member_vec.size());
        cluster_pip.push_back(cluster_prob);
        if(cluster_prob >= cluster_pip_thresh){
            grp_vec.push_back(member_vec);
            grpr2_map[i] = grp_vec.size()-1;
        }
        //printf("cluster %d:  %7.3f\n\n", i+1, cluster_prob);
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

    for(int i=0;i<nsnp_vec_sort.size();i++){
        if(nsnp_vec_sort[i].incl_prob<1e-3&&!output_all)
            break;
        if(nsnp_vec_sort[i].cluster==-1&&!output_all)
            continue;
        fprintf(outfd,"((%d))\t%15s\t%8.5e\t%7.3f\t%d\n",i+1, nsnp_vec_sort[i].name.c_str(), nsnp_vec_sort[i].incl_prob, single_log10_abfv[nsnp_vec_sort[i].name], nsnp_vec_sort[i].cluster);
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

}



void controller::parse_nmodel(Nmodel nmod){

    istringstream iss(nmod.id);
    string token;
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
        nsnp_vec[index].incl_prob += nmod.prob;

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

    for(int i=0;i<p;i++){

        vector<int>  mcfg = null_config;
        mcfg[i] = 1;
        fprintf(outfd,"%25s %10s   %9.3f\n" , pars.geno_map[i].c_str(), pars.pheno_name.c_str(), mlr.compute_log10_ABF(mcfg));

    }  

}


void controller::extract_ss(){

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
        fprintf(fd,"%15s\t%7.3e\n",name.c_str(), gsl_matrix_get(Z,i,0));
    }
    fclose(fd);

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
    if(!use_ss){
        double *gi = &pars.geno_vec[0][i][0];
        double *gj = &pars.geno_vec[0][j][0];
        r2 = pow(gsl_stats_correlation(gi,1,gj, 1, pars.geno_vec[0][i].size()),2);
    }else{
        r2 = pow(gsl_matrix_get(pars.ld_matrix,i,j),2);
    }
    return r2;
}




void size_model::update(){

    map<string, double>::iterator iter;

    vector<double> postv;
    for(iter = post_map.begin(); iter != post_map.end(); iter++){
        postv.push_back(iter->second);
    }

    vector<double> wv(postv.size(),1.0);
    log10_sum_post = log10_weighted_sum(postv,wv);

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


