using namespace std;

#include "parser.h"
#include "stdlib.h"
#include <string>
#include <sstream>
#include <fstream>

void parser::process_data(char *filename){
  
  
  // reading data file
  ifstream infile(filename);
  string line;

  while(getline(infile,line))
    process_line(line);
  
  //output();
    
}


void parser::process_line(string line){
    
  // parsing

  istringstream ins(line);
  
  char *content = new char[strlen(line.c_str())+1];
  strcpy(content, line.c_str());


  string header;
  ins>>header;

  string name;
  string grp;
  ins>> name;
  ins >> grp;
  

  vector<double> vecv;
  
  int missing_flag = 0;  // missing data flag
  string val;
  while(ins>>val){
    
    if(val == "NA" || val == "na"){
      missing_flag = 1;
      vecv.push_back(-1);
    }else
      vecv.push_back(atof(val.c_str()));
  }
  
  // impute missing data with mean

  if(missing_flag == 1){
    
    double nm_sum = 0;
    int nm_count = 0;
    vector<int> iv;
    for (int i=0;i<vecv.size();i++){
      if(vecv[i]==-1){
	iv.push_back(i);
      }else{
	nm_sum += vecv[i];
	nm_count++;
      }
    }

    double imp_mean = 0;
    if(nm_count!=0){
      imp_mean = nm_sum/nm_count;
    }

    for(int i=0;i<iv.size();i++){
      vecv[iv[i]] = imp_mean;
    }
  }
      


  if( header == "pheno" || header == "response"){
    pheno_name = name;
    pheno_vec.push_back(vecv);  
    pheno_map[pheno_vec.size()-1] = grp;
    
    pheno_index[string(grp)] = pheno_vec.size()-1;
    
    vector<vector<double> > gvec;
    geno_vec.push_back(gvec);
    vector<vector<double> > cvec;
    covar_vec.push_back(cvec);
  }else if( header == "geno" || header == "covariate" ){
    int index = pheno_index[grp];
    geno_vec[index].push_back(vecv);
    if(geno_rmap.find(name)==geno_rmap.end()){
      geno_map[geno_vec[index].size()-1] = name;
      geno_rmap[name] = geno_vec[index].size()-1;
    }
  }else if(header == "controlled"){
    int index = pheno_index[grp];
    covar_vec[index].push_back(vecv);
    
  }
  else{
    fprintf(stderr, "unknown header %s ... skipping \n", header.c_str());
  }




}

  
void parser::output(){
  
  for(int i=0;i<pheno_vec.size();i++)
    printf("%10s  ", pheno_map[i].c_str());
    
  
  for(int i=0;i<geno_vec[0].size();i++)
    printf("%10s  ", geno_map[i].c_str());

  printf("\n");

  int n = pheno_vec[0].size();
  
  for(int i=0;i<n;i++){
    for(int j=0;j<pheno_vec.size();j++)
      printf("%9.3f  ",pheno_vec[j][i]);
    
    for(int j=0;j<geno_vec.size();j++){
      printf("  ");
      for(int k=0;k<geno_vec[j].size();k++){
	printf("%.1f ", geno_vec[j][k][i]);
      }
      printf("\n");
    }
  }
}
