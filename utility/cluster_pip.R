args = commandArgs(trailingOnly = TRUE)
file = args[1]
cut = 0

if(length(args)==2){
  cut = as.numeric(args[2])
}

d = read.table(file,head=T)
attach(d)
dm = as.dist(d[,2:dim(d)[2]])
cm = hclust(dm)


if(cut==0){

  ms = mean(sum(prob))

  if(ms<=1){
    gp1 = cutree(cm,1)
    v = sum(prob)
    cluster = cbind(prob,gp1)
  }


  if(ms>1 & ms<2){
    gp1 = cutree(cm,1)
    gp2 = cutree(cm,2)

    cluster1 = cbind(prob,gp1)
    cluster2 = cbind(prob,gp2)

    v1 = sum(prob)
    v2 = sapply(1:2,function(x) sum(prob[gp2==x]))
  
    if(v1<=1.01){
      v= v1
      cluster = cluster1
    }else{
      v= v2
      cluster = cluster2
   }
 

  }

 
  if(ms>=2){
    K1 = floor(sum(prob))
    K2 = ceiling(sum(prob))

    gp1 = cutree(cm,K1)
    gp2 = cutree(cm,K2)

    cluster1 = cbind(prob,gp1)
    cluster2 = cbind(prob,gp2)

    v1 = sapply(1:K1,function(x) sum(prob[gp1==x]))

    v2 = sapply(1:K2,function(x) sum(prob[gp2==x]))

    if(var(v1)<=var(v2)){
      v= v1
      cluster = cluster1
    }else{
      v= v2
      cluster = cluster2
   }
 }

 


}

if(cut>0){

   gp1 = cutree(cm,cut)
   v = sapply(1:cut, function(x) sum(prob[gp1==x]))
   cluster = cbind(prob,gp1)
}




probv = cbind(1:length(v), v)
#colnames (probv) = c("Cluster", "PIP")
#rownames (probv) = rep(" ",length(v))
cluster= cluster[order(cluster[,2]),]
probc = sapply(cluster[,2], function(x) probv[x,2])
cluster = cbind(cluster, probc)
colnames(cluster) = rep("",3)
cluster
#probv

