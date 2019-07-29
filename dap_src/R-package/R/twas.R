myRegression = function(formula, data, weight){
  my.lm = lm(formula, data)
  my.coef = weight * coef(my.lm)[-1]
  er2 = weight * summary(my.lm)$r.squared
  return(list(coef = my.coef, er2 = er2))
}
