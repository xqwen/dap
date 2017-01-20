using namespace std;

#include "controllerZ.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
  

  char data_file[128];
  //char ld_file[128];
  char prior_file[128];  
  char subset_file[128];
  
  //memset(grid_file,0,128);
  memset(data_file,0,128);
  memset(subset_file,0, 128);
  memset(prior_file,0,128);

  int load_bf = 0;

  
  double pi1 = -1;
  
  
  for(int i=1;i<argc;i++){
     

    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      strcpy(data_file,argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "--load_bf")==0){
      load_bf = 1;
      continue;
    }


    /*
    if(strcmp(argv[i], "-ld")==0 || strcmp(argv[i], "-LD")==0){
      strcpy(ld_file,argv[++i]);
      continue;
    }
    */
    
    
    if(strcmp(argv[i], "-subset")==0){
      strcpy(subset_file,argv[++i]);
      continue;
    }


    // prior file
    if(strcmp(argv[i], "-prior")==0 || strcmp(argv[i], "-p")==0 ){
      strcpy(prior_file,argv[++i]);
      continue;
    }

    
    // prior options

    if(strcmp(argv[i], "-pi1")==0 ){
      pi1 = atof(argv[++i]);
      continue;
    }
        
    fprintf(stderr, "Error: unknown option \"%s\"\n",argv[i]);
    exit(1);
   
  }
  
  controllerZ con;

  if(strlen(subset_file)!=0){
    con.set_subset(subset_file);
  }

  if(load_bf==1)
    con.init_BF(data_file);
  else
    con.init(data_file);
  
  
  if(strlen(prior_file)==0)
    con.set_prior(pi1);
  else
    con.set_prior(prior_file);
  
  con.run();
  return 1;
   

}


