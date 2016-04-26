# Fit and returns the climate response decompositon that will be used as the covariate x
# E[x] = beta_nat * nat + beta_ghg * ghg + beta_others * others
boot_lm_gno <- function(psamples){
  l_lm <- lapply(psamples$bsamples, lm_gno)
  l_lm_an  <- mapply(predict_gno, l_lm, newdata=psamples$bsamples, SIMPLIFY=FALSE)
  l_lm_an_origin  <- mapply(predict_gno, l_lm,  newdata=psamples$osamples, SIMPLIFY=FALSE)
  list("l_lm"=l_lm, "l_x_an"=l_lm_an, "l_x_an_origin"=l_lm_an_origin)
}

# Fit on one time series
lm_gno <- function(data){
  lm_fit <- lm(y ~ ghg + nat + others, data=data)
  lm_fit
}


