#include "classdef.h"
#include <fstream>
#include <sstream>



void show_banner(){

  fprintf(stderr, "\n\nTORUS: QTL Discovery Utilizing Genomic Annotations\n\n");
  fprintf(stderr, "Usage: torus -est|qtl|dump_prior  -d input_data.gz [-smap snp_map.gz] [-gmap gene_map.gz] [-annot annotation_file.gz] [--load_bf | --load_zval]\n\n");
  


}



int main(int argc, char **argv){
  
  // creating the grid
  
  //olist.push_back(0.1);
  //phlist.push_back(0.05);

  char data_file[128];
  char gmap_file[128];
  char smap_file[128];
  char annot_file[128];
  char prior_dir[128];
  char output_pip[128];

  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  
  memset(data_file,0,128);


  char init_file[128];


  int find_egene = 0;
  int est = 0;
  double alpha = 0.05;

  double init_pi1 = 1e-3;

  char ci_file[128];
  memset(ci_file,0,128); 
  memset(data_file,0,128); 
  memset(gmap_file,0,128);
  memset(smap_file,0,128);
  memset(annot_file,0,128);
  memset(init_file,0,128);
  memset(prior_dir,0,128);
  memset(output_pip,0,128);

  int force_logistic = 0;
  int prob_annot = 0;
  
  int data_format = 1;

  double EM_thresh = 0.05;
  double dist_bin_size = -1;
  double l1_lambda = 0;
  double l2_lambda = 0;
  
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
    
    if(strcmp(argv[i], "-l1_lambda")==0){
      l1_lambda = atof(argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-l2_lambda")==0){
      l2_lambda = atof(argv[++i]);
      continue;
    }
    

    if(strcmp(argv[i], "--single_fuzzy_annot")==0){
      prob_annot = 1;
      continue;
    }

    
    if(strcmp(argv[i], "-init_pi1")==0){
      init_pi1 = atof(argv[++i]);
      continue;
    }


    if(strcmp(argv[i], "--force_logistic")==0){
      force_logistic = 1;
      continue;
    }
    
    if(strcmp(argv[i], "--load_bf")==0 || strcmp(argv[i], "--bf")==0 ){
      data_format = 2;
      continue;
    }
    
   
    if(strcmp(argv[i], "--load_zval")==0 || strcmp(argv[i], "--zval")==0){
      data_format = 3;
      continue;
    }
       
    
    if(strcmp(argv[i], "--load_fastqtl")==0 || strcmp(argv[i], "--fastqtl")==0){
      data_format = 4;
      continue;
    }

    if(strcmp(argv[i], "--load_matrixeqtl")==0 || strcmp(argv[i], "--matrixeqtl")==0){
      data_format = 1;
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
    
    
    if(strcmp(argv[i], "-egene")==0 || strcmp(argv[i], "-qtl")==0 ){
      find_egene = 1;
      continue;
    }

    if(strcmp(argv[i], "-dump_prior")==0){
      strcpy(prior_dir, argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-dump_pip")==0){
      strcpy(output_pip, argv[++i]);
      continue;
    }


    if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0 ){
      show_banner();
      continue;
    }



    if(strcmp(argv[i], "-alpha")==0){
      alpha = atof(argv[++i]);
    }


    fprintf(stderr, "Error: undefined option %s\n", argv[i]);
    show_banner();
    exit(0);

  }    



  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");
    show_banner();
    exit(0);
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

  if(prob_annot){
    con.single_fuzzy_annot = 1;
  }

  if(l1_lambda!=0){
    con.l1_lambda = l1_lambda;
    con.force_logistic = 1;
  }

  if(l2_lambda!=0){
    con.l2_lambda = l2_lambda;
    con.force_logistic =1;
  }
  
  
  con.init_pi1 = init_pi1;


  switch(data_format){
  case 1:
    con.load_data(data_file);
    break;
  case 2:
    con.load_data_BF(data_file);
    break;
  case 3:
    con.load_data_zscore(data_file);
    break;
  case 4:
    con.load_data_fastqtl(data_file);
    gmap_file[0]=smap_file[0] = 0;
    break;
  default:
    con.load_data(data_file);
    break;
  }
   

  con.load_map(gmap_file, smap_file);  
  con.load_annotation(annot_file);
  fprintf(stderr,"Initializing ... \n");

  if(est==0 && find_egene==0){
    est = 1;
  }


  if(est)
    con.estimate();
  if(find_egene){
    con.find_eGene(alpha);
  }
  if(strlen(prior_dir)>0){
    con.dump_prior(prior_dir);
  }
  if(strlen(output_pip)>0){
    fprintf(stderr,"#-dump_pip output_pip file: %s \n",output_pip);
    con.dump_pip(output_pip);
  }
}


