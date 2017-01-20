using namespace std;

#include "controller.h"
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv){
  
  char grid_file[128];
  char data_file[128];
  char out_file[128];
  char gene_name[64];
  char inc_file[128];

  char prior_file[128];
  
  
  memset(gene_name,0,64);
  memset(grid_file,0,128);
  memset(data_file,0,128);
 
  memset(out_file,0,128);
  memset(inc_file,0,128);
  memset(prior_file,0,128);
  
  
  double abf_option = 0.5;
  int abf_prior_option = 1;

  int msize =-1;
  
  int set_yty = 0;
  double yty;
  
  int sample_size = -1;

  int output_all = 0;

  double pi1 = -1;
  double lambda = 0.5;

  
  double snp_select_thresh = -1;
  double size_select_thresh = -1;

  int est_option = 0;

  int thread = 1;
  
  for(int i=1;i<argc;i++){
     

    // required data files and additional info

    if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-grid")==0){
      strcpy(grid_file,argv[++i]);
      continue;
    } 
    
    if(strcmp(argv[i], "-s")==0 || strcmp(argv[i], "-summary")==0){
      strcpy(data_file,argv[++i]);
      continue;
    }


    
    if(strcmp(argv[i], "-n")==0 || strcmp(argv[i], "-N")==0){
      sample_size = atoi(argv[++i]);
      continue;
    }



    if(strcmp(argv[i],"-include")==0 || strcmp(argv[i],"-inc")==0){
      strcpy(inc_file, argv[++i]);
      continue;
    }
								   


    if(strcmp(argv[i], "-yty")==0 ){
      yty = atof(argv[++i]);
      set_yty = 1;
      continue;
    }


    // prior file
    if(strcmp(argv[i], "-prior")==0 || strcmp(argv[i], "-p")==0 ){
      strcpy(prior_file,argv[++i]);
      continue;
    }

    
    // output file
    if(strcmp(argv[i], "-o")==0){
      strcpy(out_file, argv[++i]);
      continue;
    }


    // prior options

    if(strcmp(argv[i], "-pi1")==0 ){
      pi1 = atof(argv[++i]);
      continue;
    }
    
    // abf option
    if(strcmp(argv[i],"-abf")==0){
      abf_option = atof(argv[++i]);
      continue;
    }
    
  

    // thresholds
    
    if(strcmp(argv[i],"-size_thresh")==0 || strcmp(argv[i],"-st")==0 ){
      size_select_thresh = atof(argv[++i]);
      continue;
    }

    if(strcmp(argv[i],"-inc_thresh")==0 || strcmp(argv[i],"-it")==0 ){
      snp_select_thresh = atof(argv[++i]);
      continue;
    }


   
    // msize option 
    if(strcmp(argv[i],"-msize")==0){
      msize = atoi(argv[++i]);
      continue;
    }


    // openmp threads

    if(strcmp(argv[i], "-t")==0){
      thread = atoi(argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "-all")==0){
      output_all = 1;
      continue;
    }


    
    fprintf(stderr, "Error: unknown option \"%s\"\n",argv[i]);
    exit(1);
   
  }
  
  controller con;
  

  // set options
  if(sample_size == -1){
    fprintf(stderr, "Error: sample size is not specified (-n sample_size)\n");
    exit(1);
  }


  con.set_sample_size(sample_size);

  if(set_yty == 1){
    con.set_yty(yty);
  }

  con.set_outfile(out_file);
  con.set_gene(gene_name);
  con.set_sslr_option(abf_option);
  con.set_thread(thread);

   
  if(output_all == 1)
    con.set_output_all();

  if(snp_select_thresh>=0)
    con.set_snp_select_thresh(snp_select_thresh);
  
  if(size_select_thresh >=0)
    con.set_size_select_thresh(size_select_thresh);
  
  

  // init
  
  
 
  con.init(grid_file, data_file, inc_file);
  
  
  if(strlen(prior_file)==0)
    con.set_prior(pi1,lambda);
  else
    con.set_prior(prior_file);
  
  con.init_prior();
  
  // run
  con.run();
  return 1;
   

}


