#include "classdef.h"
#include <fstream>
#include <sstream>



int main(int argc, char **argv){
  
  // creating the grid
  
  //olist.push_back(0.1);
  //phlist.push_back(0.05);

  char data_file[128];
  char gmap_file[128];
  char smap_file[128];
  char annot_file[128];
  char prior_dir[128];

  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  memset(data_file,0,128);


  char init_file[128];


  int find_egene = 0;
  int est = 1;
  

  

  char ci_file[128];
  memset(ci_file,0,128); 
  memset(data_file,0,128); 
  memset(gmap_file,0,128);
  memset(smap_file,0,128);
  memset(annot_file,0,128);
  memset(init_file,0,128);
  memset(prior_dir,0,128);
  

  int force_logistic = 0;
  int load_bf = 0;
  double EM_thresh = 0.05;
  double dist_bin_size = -1;

  for(int i=1;i<argc;i++){
    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      strcpy(data_file,argv[++i]);
      continue;
    }



    if(strcmp(argv[i], "-gmap")==0){
      strcpy(gmap_file,argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "-smap")==0 ){
      strcpy(smap_file,argv[++i]);
      continue;
    }


    if(strcmp(argv[i], "-annot")==0 ){
      strcpy(annot_file,argv[++i]);
      continue;
    }
    
    if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-thresh")==0){
      EM_thresh = atof(argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "--force_logistic")==0){
      force_logistic = 1;
      continue;
    }
    
    if(strcmp(argv[i], "--load_bf")==0){
      load_bf = 1;
      continue;
    }
    

    if(strcmp(argv[i], "-dist_bin_size") == 0){
      dist_bin_size = atof(argv[++i]);
      continue;
    }

    
    
    if(strcmp(argv[i], "-est")==0){
      est = 1;
      continue;
    }
    
    
    if(strcmp(argv[i], "-egene")==0){
      find_egene = 1;
      continue;
    }

    if(strcmp(argv[i], "-dump_prior")==0){
      strcpy(prior_dir, argv[++i]);
      continue;
    }


    //fprintf("Error: undefined option %s\n", argv[i]);

  }    



  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");
    exit(1);
  }

    
  // a global variable 
  controller con;
  con.EM_thresh = EM_thresh;
  
  if(dist_bin_size > 0){
    con.dist_bin_size = dist_bin_size;
  }

  if(force_logistic){
    con.force_logistic = 1;
  }



  if(load_bf){
    con.load_data_BF(data_file);
  }else{
    con.load_data(data_file);
  }
  con.load_map(gmap_file, smap_file);
  con.load_annotation(annot_file);
  
  fprintf(stderr,"Initializing ... \n");
  if(est)
    con.estimate();
  if(find_egene)
    con.find_eGene();
  
  if(strlen(prior_dir)>0){
    con.dump_prior(prior_dir);
  }
}


