
sim_pheno<-function(d){
  int = runif(1, min=-2,max=2)
  n = dim(d)[2]-3
  

  if(nsig==0){
    yv = int + rnorm(n)
  }else{
    G = as.matrix(d[index,4:dim(d)[2]])
    beta = apply(as.matrix(bbar),1, function(x) rnorm(1,mean=x,sd=abs(bbar)*0.1))
    #beta = matrix(ncol=1,rep(0.6,length(bbar)))
    yv = int + t(G)%*%beta + rnorm(n)	
  }
  return(yv)
} 


args = commandArgs(trailingOnly = TRUE)
gene = args[1]
lambda = as.numeric(args[2])
ceu_data = paste("geno_data/",gene,".ceu.geno",sep="")
fin_data = paste("geno_data/",gene,".fin.geno",sep="")
gbr_data = paste("geno_data/",gene,".gbr.geno",sep="")
tsi_data = paste("geno_data/",gene,".tsi.geno",sep="")
yri_data = paste("geno_data/",gene,".yri.geno",sep="")



ceu_out = paste("pheno_data/",gene,".ceu.pheno",sep="")
fin_out = paste("pheno_data/",gene,".fin.pheno",sep="")
gbr_out = paste("pheno_data/",gene,".gbr.pheno",sep="")
tsi_out = paste("pheno_data/",gene,".tsi.pheno",sep="")
yri_out = paste("pheno_data/",gene,".yri.pheno",sep="")



annot= rep(0,50)
set=1:50
annot[which(set%%5==0)]=1
logit=-4+annot*lambda
prob = exp(logit)/(1+exp(logit))
index = which(sapply(prob, function(x) rbinom(1,1,x))==1)
nsig = length(index)



d1 = read.table(ceu_data)
d2 = read.table(fin_data)
d3 = read.table(gbr_data)
d4 = read.table(tsi_data)
d5 = read.table(yri_data)

if(nsig !=0){
  snp_set = as.character(d1[index,2])
  truth_out =  paste("truth/",gene,".truth",sep="")
  bbar = rnorm(nsig,sd=0.6)
  td = cbind(snp_set,bbar)
  write(file=truth_out,t(td),ncol=2)
}

annot_out = paste("truth/",gene,".annot",sep="")
ad = cbind(as.character(d1[,2]),annot)
write(file=annot_out,t(ad),ncol=2)


yv1 = sim_pheno(d1)
write(file=ceu_out,ncol=1,yv1)


yv2 = sim_pheno(d2)
write(file=fin_out,ncol=1,yv2)
yv3 = sim_pheno(d3)
write(file=gbr_out,ncol=1,yv3)
yv4 = sim_pheno(d4)
write(file=tsi_out,ncol=1,yv4)
yv5 = sim_pheno(d5)
write(file=yri_out,ncol=1,yv5)






