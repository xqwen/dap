attach(read.table("regress.dat"))
m = glm(V4~V3,family="binomial")
summ = summary(m)
mu = summ$coef[1,1]
lambda = summ$coef[2,1]

logit = exp(mu+lambda*V3)
prob = logit/(1+logit)
dv = cbind(as.character(V1),as.character(V2),prob)


c(summ$coef[1,1], summ$coef[2,1])
write(file="prob.out",t(dv),ncol=3)


