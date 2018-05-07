tot = 4613
fdr_control<-function(lfdr,truth,alpha=0.05){
    lfdr_sort = sort(lfdr)
    fdr = cumsum(lfdr_sort)/(1:length(lfdr_sort))
    index = max(which(fdr<=alpha))
    lfdr_thresh = lfdr_sort[index]
    fd = length(truth[truth==0&lfdr<=lfdr_thresh])
    td = length(truth[truth==1&lfdr <= lfdr_thresh])
    fdp = fd/index
    power = td/tot
    return(c(fdp,power))
}
 

d1 = read.table("dap_sss.cluster_pip.dat")
lfdr1 = 1-d1$V3
truth1 = d1$V4

d2 = read.table("dap_z.cluster_pip.dat")
lfdr2 = 1-d2$V3
truth2 = d2$V4


thresh = c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25)
t(sapply(thresh,function(x) fdr_control(lfdr1,truth1,x)))

t(sapply(thresh,function(x) fdr_control(lfdr2,truth2,x)))

