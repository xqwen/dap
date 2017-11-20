attach(read.table("logit.dat"))
m = summary(glm(V3~V2,family="binomial"))
beta = m$coef[2,1]
sd = m$coef[2,2]
c(m$coef[1,1],beta, beta-1.96*sd, beta+1.96*sd)