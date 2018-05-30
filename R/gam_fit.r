#' @importFrom mgcv gam 
#' @importFrom mgcv predict.gam 
#' @importFrom mgcv s
#' @importFrom reshape2 melt 
#' @importFrom reshape2 dcast 

# from dataset with ebm simulated responses (with random scaling factors or not), create bootstrap datasets 
# and corresponding original/complete datasets
prepare_osamples <- function(osamples, bindexes){
  bsamples <- mapply(function(data, ind){ans <- data[ind, ]
                     ans[order(ans$time), ]},
                     osamples, bindexes, SIMPLIFY=FALSE)
  osamples <- lapply(osamples, compute_original_data)
  list("bsamples"=bsamples, "osamples"=osamples)
}

# Fit and returns the climate response decompositon that will be used as the covariate x
# E[x] = beta_nat * nat + s(ant)
boot_gam_allnat <- function(psamples, knots=NULL, fixed=FALSE){
  # model fit
  l_gam <- lapply(psamples$bsamples, function(data) gam_allnat(data, fixed=fixed, knots=knots))
  # predict responses on bootstrap samples
  l_gam_an  <- mapply(predict_allnat, l_gam, newdata=psamples$bsamples, SIMPLIFY=FALSE)
  # predict responses on original/complete data
  l_gam_an_origin  <- mapply(predict_allnat, l_gam,  newdata=psamples$osamples, SIMPLIFY=FALSE)
  list("l_gam"=l_gam, "l_x_an"=l_gam_an, "l_x_an_origin"=l_gam_an_origin)
}

# Fit and returns the climate response decompositon that will be used as the covariate x
# E[x] = beta_nat * nat + beta_ghg * ghg + s(others)
boot_gam_gno <- function(psamples, knots=NULL, fixed=FALSE){
  # model fit
  l_gam <- lapply(psamples$bsamples, function(data) gam_gno(data, fixed=fixed, knots=knots))
  # predict responses on bootstrap samples
  l_gam_an  <- mapply(predict_gno, l_gam, newdata=psamples$bsamples, SIMPLIFY=FALSE)
  # predict responses on original/complete data
  l_gam_an_origin  <- mapply(predict_gno, l_gam,  newdata=psamples$osamples, SIMPLIFY=FALSE)
  list("l_gam"=l_gam, "l_x_an"=l_gam_an, "l_x_an_origin"=l_gam_an_origin)
}

# Fit on one time series
gam_allnat <- function(data, knots=NULL, fixed=FALSE ){
  if(!is.null(knots))
    # gam <- fit <- gam(x ~ s(ant, k=knots, fx=fixed, bs="ad",
    # xt=list(bs="cs")) + nat, data=data, control=list(keepData=TRUE))
    gam_fit <- gam(x ~ s(ant, k=knots, fx=fixed, bs="ts") + nat, data=data)#, control=list(keepData=TRUE))
  else
    gam_fit <- gam(x ~ s(ant, bs="ts") + nat, data=data)#, control=list(keepData=TRUE))
  gam_fit
}

# Fit on one time series
gam_gno <- function(data, knots=NULL, fixed=FALSE ){
  if(!is.null(knots))
    # gam <- fit <- gam(x ~ s(ant, k=knots, fx=fixed, bs="ad",
    # xt=list(bs="cs")) + nat, data=data, control=list(keepData=TRUE))
    gam_fit <- gam(x ~ s(others, k=knots, fx=fixed, bs="ts") + nat + ghg, data=data)#, control=list(keepData=TRUE))
  else
    gam_fit <- gam(x ~ s(others, bs="ts") + nat + ghg, data=data)#, control=list(keepData=TRUE))
  gam_fit
}

# Given a fit, returns the three components: x_all=x_nat+x_ant
predict_allnat <- function(gam_fit, newdata){
  gam_all <- predict(gam_fit, newdata=newdata, type="response")
  data_nat <- newdata
  data_nat$ant <- 0
  data_ant <- newdata
  data_ant$nat <- 0
  gam_nat <- predict(gam_fit, newdata=data_nat, type="response")
  gam_ant <- predict(gam_fit, newdata=data_ant, type="response")
  cbind(newdata, "x_all"=gam_all, "x_nat"=gam_nat, "x_ant"=gam_ant)
}

# Given a fit, returns the three components: x_all=x_nat+x_ant
predict_gno <- function(gam_fit, newdata){
  gam_all <- predict(gam_fit, newdata=newdata, type="response")
  data_nat <- newdata
  data_nat$ghg <- 0
  data_nat$others <- 0
  data_ant <- newdata
  data_ant$nat <- 0
  gam_nat <- predict(gam_fit, newdata=data_nat, type="response")
  gam_ant <- predict(gam_fit, newdata=data_ant, type="response")
  cbind(newdata, "x_all"=gam_all, "x_nat"=gam_nat, "x_ant"=gam_ant)
}

# takes a dataset, remove x and y, remove historical nat runs and doublons, and order by time
compute_original_data <- function(original_data){
  names_od <- names(original_data)
  original_data <- original_data[, !(names_od %in% c("y", "x"))]
  # original_data <- subset(original_data, hnat==0)
  original_data <- unique(original_data)
  original_data <- original_data[order(original_data$time), ]
  original_data
}

# compute means of ant covariates in the pre-industrial
compute_ant_bias <- function(bsamples, ant_var, pre_ind=c(1850, 1879)){
    bsamples <- lapply(bsamples, function(data) subset(data, time >= pre_ind[1] & time <= pre_ind[2], select=ant_var))
    bbias <- lapply(bsamples, function(data) apply(data, 2, mean))
    return(bbias)
}

# shift ANT covariates so that their means in the pre-industrial is zero
# Can either compute the mean for each bootstrap sample or only one the first one 
boot_correct_ant_bias <- function(psamples, ant_var,  pre_ind=c(1850, 1879), reuse_ant_bias=FALSE){
  bsamples <- psamples$bsamples
  osamples <- psamples$osamples
  if(reuse_ant_bias){
    bbias <- unlist(compute_ant_bias(bsamples[1], ant_var=ant_var, pre_ind=pre_ind))
    correct_ant_bias_reuse <- function(data){
      data[, ant_var]  <- sweep(data[, ant_var, drop=FALSE], 2, bbias)
      data
    }
    bsamples <- lapply(bsamples, correct_ant_bias_reuse)
    osamples <- lapply(osamples, correct_ant_bias_reuse)
  } else{
    bbias <- compute_ant_bias(bsamples, ant_var=ant_var, pre_ind=pre_ind)
    correct_ant_bias <- function(data, bias){
      data[, ant_var]  <- sweep(data[, ant_var, drop=FALSE], 2, bias)
      data
    }
    bsamples <- mapply(correct_ant_bias, data=bsamples, bias=bbias, SIMPLIFY=FALSE)
    osamples <- mapply(correct_ant_bias, data=osamples, bias=bbias, SIMPLIFY=FALSE)
  }
  list("bsamples"=bsamples, "osamples"=osamples)
}


