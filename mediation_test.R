## a small expriment to test independence and correlation

s = rnorm(1000, mean = 10, sd = 1)
expr = rnorm(1000, mean = 10, sd = 1)
error = function() {rnorm(1000, mean = 0, sd = 1)}

#model 1: smoking ---> methylation <-+- expression
methy = s*1.5 - expr*1.5 + error()

summary(lm(expr~ s))
summary(lm(expr~ s + methy))
summary(lm(expr ~ methy))
summary(lm(methy~ s))
summary(lm(methy~ s + expr))

## model 2: smoking ---> methylation -+->expression ; smoking -+-> expression

methy = s*(-0.2) + error()
expr = s*0.5 + methy*(0.3) + error()

summary(lm(expr ~ s))
summary(lm(methy ~ s + expr))
summary(lm(expr ~ methy))
summary(lm(methy ~ s))
