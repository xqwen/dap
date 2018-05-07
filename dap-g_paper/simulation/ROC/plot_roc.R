pdf(file="roc.pdf",width=8,height=8,bg="white")

d1 = read.table("dap_sss.roc.sort.dat")
tv= sort(unique(d1$V2),decreasing=T)
ttv = tv[tv>1e-4]
rst1  = t(sapply(ttv, function(x) c(sum(d1$V3[d1$V2>=x]), length(d1$V3[d1$V2>=x]))))
fd1 = rst1[,2]-rst1[,1]
plot(rst1[,1]~fd1,t="l",col="brown",xlim=c(0,3000),ylab = "True discoveries", xlab ="False discoveries")


d2 = read.table("dap_z.roc.sort.dat")
tv= sort(unique(d2$V2),decreasing=T)
ttv = tv[tv>1e-4]
rst2  = t(sapply(ttv, function(x) c(sum(d2$V3[d2$V2>=x]), length(d2$V3[d2$V2>=x]))))
fd2 = rst2[,2]-rst2[,1]
lines(rst2[,1]~fd2,col="darkgreen")


d3 = read.table("finemap.roc.sort.dat")
tv= sort(unique(d3$V2),decreasing=T)
ttv = tv[tv>1e-4]
rst3  = t(sapply(ttv, function(x) c(sum(d3$V3[d3$V2>=x]), length(d3$V3[d3$V2>=x]))))
fd3 = rst3[,2]-rst3[,1]
lines(rst3[,1]~fd3,col="navyblue")



d4 = read.table("single.roc.sort.dat")
tv= sort(unique(d4$V2),decreasing=T)
ttv = tv[tv>-0.1]
rst4  = t(sapply(ttv, function(x) c(sum(d4$V3[d4$V2>=x]), length(d4$V3[d4$V2>=x]))))
fd4 = rst4[,2]-rst4[,1]
lines(rst4[,1]~fd4,col="magenta")


legend("bottomright", legend=c("DAP-G (sufficient summary stats)", "DAG-G (z-scores)", "FINENAMP (z-scores)", "Single SNP Testing"),
       col=c("brown", "darkgreen", "navyblue", "magenta"), lty=c(1,1,1,1), cex=0.8)

dev.off()

