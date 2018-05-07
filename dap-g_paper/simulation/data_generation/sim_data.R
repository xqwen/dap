args = commandArgs(trailingOnly=TRUE)
index = args[1]

X = as.matrix(read.table("sim.r.geno"))
n = dim(X)[1]      # number of observations
p = dim(X)[2]      # number of variables
freq = 0.05

region = rbinom(91,1,freq)
amplitude = 0.6

get.beta<-function(v){
  rstv = rep(0,11)
  if(v==1){
    index = sample(1:11,1)
    rstv[index]=rnorm(1,sd=amplitude);

  }
  return(rstv)
}

beta = as.vector(sapply(region, get.beta))

y.sample <- function() X %*% beta + rnorm(n)

Y = matrix(ncol=1,y.sample())
dY = c("pheno gene group", as.vector(Y))
rsv = paste("geno rs",seq(1,p),"  group", sep="")
dX = cbind(rsv, t(X))

data = rbind(dY,dX)
filename = paste("sim.",index,".dat", sep="")
write(file=filename, t(data), ncol = dim(data)[2])

rdata = cbind(Y,X)
filename = paste("sim.",index,".r.dat", sep="")
write(file=filename,t(rdata), ncol=dim(rdata)[2])

tdata = cbind(rsv, beta)[beta!=0,]
filename = paste("sim.",index,".truth", sep="")
write(file=filename,t(tdata),ncol=2)

