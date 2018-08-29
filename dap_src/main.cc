using namespace std;

#include "controller.h"
#include <stdlib.h>

int main(int argc, char **argv){

    char grid_file[128];
    char data_file[128];
    char zval_file[128];
    char est_file[128];
    char ld_file[128];
    char out_file[128];
    char log_file[128];
    char gene_name[64];

    char prior_file[128];


    memset(gene_name,0,64);
    memset(grid_file,0,128);
    memset(out_file,0,128);
    memset(log_file,0,128);
    memset(data_file,0,128);
    memset(ld_file,0,128);
    memset(zval_file,0,128);
    memset(est_file,0,128);

    memset(prior_file,0,128);


    int ld_format = 1; // for correlation matrix


    double abf_option = 0.5;
    int abf_prior_option = 1;

    int msize =-1;

    int output_all = 0;

    double pes = 1.0;
    double pi1 = -1;
    double lambda = 0.5;
    double ld_control = -1;

    int size_limit = -1;

    int sample_size = -1;
    double syy = -1;

    double snp_select_thresh = -1;
    double size_select_thresh = -1;

    // alternative non-fm running options
    int run_scan = 0;
    int extract_ss = 0;

    int thread = 1;



    for(int i=1;i<argc;i++){


        // required data files and additional info

        if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-grid")==0){
            strcpy(grid_file,argv[++i]);
            continue;
        } 

        if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
            strcpy(data_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-d_z")==0 || strcmp(argv[i], "-d_zval")==0 || strcmp(argv[i], "-data_zval")==0){
            strcpy(zval_file,argv[++i]);
            continue;
        } 

        if(strcmp(argv[i], "-d_ld")==0 || strcmp(argv[i], "-data_ld")==0){
            strcpy(ld_file,argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-d_ld2")==0 || strcmp(argv[i], "-data_ld")==0){
            strcpy(ld_file,argv[++i]);
            ld_format = 2;
            continue;
        }

        if(strcmp(argv[i], "-d_est")==0){ 
            strcpy(est_file,argv[++i]);
            continue;
        } 

        if(strcmp(argv[i], "-d_n")==0){
            sample_size = atoi(argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-d_syy")==0){
            syy = atof(argv[++i]);
            continue;
        }


        // prior file
        if(strcmp(argv[i], "-prior")==0 || strcmp(argv[i], "-p")==0 ){
            strcpy(prior_file,argv[++i]);
            continue;
        }


        // output file
        if(strcmp(argv[i], "-o")==0 || strcmp(argv[i], "-output")==0){
            strcpy(out_file, argv[++i]);
            continue;
        }

        if(strcmp(argv[i], "-l")==0 || strcmp(argv[i], "-logfile")==0){
            strcpy(log_file, argv[++i]);
            continue;
        }



        // prior options

        if(strcmp(argv[i], "-ens")==0 ){
            pes = atof(argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-pi1")==0){
            pi1 = atof(argv[++i]);
            continue;
        }

        // abf option
        if(strcmp(argv[i],"-abf")==0){
            abf_option = atof(argv[++i]);
            continue;
        }


        // thresholds

        if(strcmp(argv[i],"-converg_thresh")==0){
            size_select_thresh = atof(argv[++i]);
            continue;
        }

        if(strcmp(argv[i],"-size_limit")==0){
            size_limit = atoi(argv[++i]);
            continue;
        }

        if(strcmp(argv[i],"--no_size_limit")==0){
            size_limit = -1;
            continue;
        }


        if(strcmp(argv[i],"-ld_control")==0){
            ld_control = atof(argv[++i]);
            continue;
        }



        // msize option  for DAP-K
        if(strcmp(argv[i],"-msize")==0 || strcmp(argv[i], "-dapk") ==0){
            msize = atoi(argv[++i]);
            continue;
        }


        // openmp threads

        if(strcmp(argv[i], "-t")==0){
            thread = atoi(argv[++i]);
            continue;
        }


        // output option
        if(strcmp(argv[i], "--all")==0 || strcmp(argv[i],"--output_all")==0){
            output_all = 1;
            continue;
        }



        // gene/locus name
        if(strcmp(argv[i], "-name")==0 || strcmp(argv[i], "-n")==0){
            strcpy(gene_name, argv[++i]);
            continue;
        }

        // no finemappin, just sincle SNP analysis
        if(strcmp(argv[i], "--scan")==0){
            run_scan = 1;
            continue;
        }

        if(strcmp(argv[i], "--dump_summary")==0){
            extract_ss = 1;
            continue;
        }

        if(strcmp(argv[i], "--dump_summary2")==0){
            extract_ss = 2;
            continue;
        }




        fprintf(stderr, "Error: unknown option \"%s\"\n",argv[i]);
        exit(1);

    }

    controller con;
    if(strlen(data_file)!=0){
        con.initialize(data_file,grid_file);
    }else if(strlen(zval_file)!=0 && strlen(ld_file)!=0 ){
        con.initialize(zval_file, ld_file, grid_file,sample_size, ld_format);
    }else if(strlen(ld_file)!=0 && strlen(est_file)!=0 && sample_size >0 && syy>0){
        con.initialize(est_file, ld_file, grid_file, sample_size, syy,ld_format);
    }else{
        fprintf(stderr, "Error: no suitable input data specified \n");
        exit(1);
    }

    con.set_outfile(out_file, log_file);

    con.set_gene(gene_name);
    con.set_abf_option(abf_option);
    con.set_thread(thread);

    con.set_size_limit(size_limit);

    if(ld_control>=0)
        con.set_ld_control(ld_control);

    if(msize>=1){
        con.set_max_size(msize);
    }


    if(strlen(prior_file)==0){
        if(pi1 != -1){
            if(0<pi1 && pi1 <1){
                con.set_prior(pi1);
            }else{
                fprintf(stderr, "Warning: pi1 specification is outside the range, ignored...\n");
                con.set_prior_exp(pes);
            }
        }else{
            // default
            con.set_prior_exp(pes);
        }
    }else
        con.set_prior(prior_file);

    if(output_all == 1)
        con.set_output_all();

    if(snp_select_thresh>=0)
        con.set_snp_select_thresh(snp_select_thresh);

    if(size_select_thresh >=0)
        con.set_size_select_thresh(size_select_thresh);


    con.run_option = 0;

    if(run_scan){
        con.run_option = 1;
    }

    if(extract_ss==1){
        con.run_option =2;
    }
    if(extract_ss==2){
        con.run_option =3;
    }


    // all done, print all configs
    con.print_dap_config();

    con.run();
    return 1;


}


