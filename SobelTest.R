## sobel test with factor of two levels
sobel.lm = function(pred, med, out, covariates){
  NEWDAT <- data.frame(pred = pred, med = med, out = out, covariates)
  NEWDAT <- na.exclude(NEWDAT)
  covar = names(covariates)
  colnames(NEWDAT)[1:3] = c("pred", "med", "out")
  model1 <- lm(out ~ . - med, data = NEWDAT)
  model2 <- lm(out ~ ., data = NEWDAT)
  model3 <- lm(med ~ . - out, data = NEWDAT)
  mod1.out <- summary(model1)$coef
  mod2.out <- summary(model2)$coef
  mod3.out <- summary(model3)$coef
  indir <- mod3.out[2, 1] * mod2.out[3, 1]
  effvar <- (mod3.out[2, 1])^2 * (mod2.out[3, 2])^2 + (mod2.out[3, 1])^2 * (mod3.out[2, 2])^2
  serr <- sqrt(effvar)
  zvalue = indir/serr
  out <- list(`Mod1: Y~X` = mod1.out, `Mod2: Y~X+M` = mod2.out, 
              `Mod3: M~X` = mod3.out, Indirect.Effect = indir, SE = serr, 
              z.value = zvalue, N = nrow(NEWDAT))
  return(out)
}

## sobel test with factor of three levels, in which the third level is of interest
sobel.lm2 = function(pred, med, out, covariates){
  NEWDAT <- data.frame(pred = pred, med = med, out = out, covariates)
  NEWDAT <- na.exclude(NEWDAT)
  covar = names(covariates)
  colnames(NEWDAT)[1:3] = c("pred", "med", "out")
  model1 <- lm(out ~ . - med, data = NEWDAT)
  model2 <- lm(out ~ ., data = NEWDAT)
  model3 <- lm(med ~ . - out, data = NEWDAT)
  mod1.out <- summary(model1)$coef
  mod2.out <- summary(model2)$coef
  mod3.out <- summary(model3)$coef
  indir <- mod3.out[3, 1] * mod2.out[4, 1]
  effvar <- (mod3.out[3, 1])^2 * (mod2.out[4, 2])^2 + (mod2.out[4, 1])^2 * (mod3.out[3, 2])^2
  serr <- sqrt(effvar)
  zvalue = indir/serr
  out <- list(`Mod1: Y~X` = mod1.out, `Mod2: Y~X+M` = mod2.out, 
              `Mod3: M~X` = mod3.out, Indirect.Effect = indir, SE = serr, 
              z.value = zvalue, N = nrow(NEWDAT))
  return(out)
}
