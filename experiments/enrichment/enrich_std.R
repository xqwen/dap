attach(read.table("regress_std.dat"))

V4[is.na(V4)]=1
library(qvalue)
rst = qvalue(V4,fdr.level=0.001)
rej = rst$sig+0


m = glm(rej~V3,family="binomial")
summ = summary(m)
mu = summ$coef[1,1]
lambda = summ$coef[2,1]
sd = summ$coef[2,2]
c(mu, lambda, lambda-1.96*sd, lambda+1.96*sd)


