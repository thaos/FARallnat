#' @import FARg
#' @import magrittr

# Fit a statistical model using as a covariate the climate
# response to all forcings.
# 
# \code{gev_fit} fit a statistical model to relate the evolution of the
# distribution of the y-variable with respect to an estimated climate response
# of the climate to all forcings for all bootstrap samples
#
# Possible model includes a GEV, GPD and a Gaussian distribution. Functions,
# gev_fit, gpd_fit and gauss_fit from the package FARg are used to fit the
# statistical model.
# @param l_x_an boostrap samples containing the y value to be fitted, the date
# and the corresponding climate response to all forcing x_all
# @param the function used for the fit: either gauss_fit, gev_fit, gpd_fit from
# the FARg package
# @param ... additional parameters for the fit if neccessary
# @return returns a list of the fits made for all bootstrap samples
# @examples
# # Needs to be done !
fit_and_boot_allnat <- function(l_x_an, fit, ...){
  fit0 <- fit(y, l_x_an[[1]], ~x_all, ~x_all, time_var="time", ...)
  l_fit <- lapply(l_x_an, function(data) fit(y, data, ~x_all, ~x_all, time_var="time", init=fit0$par, ...))
  l_fit
}

# From a fit compute the probability of exceending xp  in an ALL, an ANT or
# NAT only world at at times l_times given the dataset original_data
get_p_allnat_onperiod <- function(object, xp, l_time, original_data,
                                  ant_nat_or_all="all"){ if(ant_nat_or_all ==
                                  "nat")
    original_data$x_all <- original_data$x_nat
  if(ant_nat_or_all == "ant")
    original_data$x_all <- original_data$x_ant
  get_p_allnat <- function(time){
    pnt <- set_pnt(time, xp, time_var=object$time_var, data=original_data)
    if(class(object) == "gpd_fit")
      get_p(object, pnt, under_threshold=TRUE)
    else
      get_p(object, pnt) 
  }
  # debug(get_p_allnat)
  print(system.time(
  ans <- sapply(l_time, get_p_allnat)
  ))
  colnames(ans) <- l_time
  ans
}

# From a fit compute the quantile corresponding to the probability of exceedance
# p in an ALL, an ANT or NAT only world at at times l_times given the dataset original_data
get_q_allnat_onperiod <- function(object, p, l_time, original_data, ant_nat_or_all="all"){
  if(ant_nat_or_all == "nat")
    original_data$x_all <- original_data$x_nat
  if(ant_nat_or_all == "ant")
    original_data$x_all <- original_data$x_ant
  get_q_allnat <- function(time){
    pnt <- set_pnt(time, p, time_var=object$time_var, data=original_data)
    get_q(object, pnt) 
  }
  print(system.time(
  ans <- sapply(l_time, get_q_allnat)
  ))
  colnames(ans) <- l_time
  ans
}


# From a fit compute the probabilites of exceending xp  in an ALL, an ANT
# and NAT only worlds at at times l_times given the dataset original_data
get_allp_onperiod <- function(object, xp, l_time, original_data){
  p_all <- get_p_allnat_onperiod(object, xp, l_time, original_data)
  p_ant <- get_p_allnat_onperiod(object, xp, l_time, original_data, ant_nat_or_all="ant")
  p_nat <- get_p_allnat_onperiod(object, xp, l_time, original_data, ant_nat_or_all="nat")
  ans <- rbind(p_all, p_ant, p_nat)
  rownames(ans) <- c(paste(rownames(p_all), "_all", sep=""), paste(rownames(p_ant), "_ant", sep=""), paste(rownames(p_nat),"_nat", sep=""))
  ans
}

# From a fit compute the quantile corresponding to the probability of exceeding
# p in an ALL, an ANT and NAT only worlds at at times l_times given the dataset original_data
get_allq_onperiod <- function(object, p, l_time, original_data){
  q_all <- get_q_allnat_onperiod(object, p, l_time, original_data)
  q_ant <- get_q_allnat_onperiod(object, p, l_time, original_data, ant_nat_or_all="ant")
  q_nat <- get_q_allnat_onperiod(object, p, l_time, original_data, ant_nat_or_all="nat")
  ans <- rbind(q_all, q_ant, q_nat)
  rownames(ans) <- c(paste(rownames(q_all), "_all", sep=""), paste(rownames(q_ant), "_ant", sep=""), paste(rownames(q_nat),"_nat", sep=""))
  ans
}

# Compute all the probabilities of exceedance for all bootstrap samples
boot_allp_onperiod  <- function(l_fit, l_x_an_origin, l_time, xp){
  l_allp <- mapply(get_allp_onperiod, object=l_fit, original_data=l_x_an_origin, MoreArgs=list("xp"=xp, "l_time"=l_time), SIMPLIFY=FALSE)
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("x_all", "x_nat", "x_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_allp2 <- mapply(rbinding, l_allp, l_x_an_origin, SIMPLIFY=FALSE)
  a_allp <- array(unlist(l_allp2), dim=c(dim(l_allp2[[1]]), length(l_allp)))
  dimnames(a_allp)[1:2] <- dimnames(l_allp2[[1]])
  dimnames(a_allp)[[3]] <- 1:length(l_allp)
  a_allp
}

# Compute all the quantileof exceeding probability p  for all bootstrap samples
boot_allq_onperiod  <- function(l_fit, l_x_an_origin, l_time, p){
  l_allq <- mapply(get_allq_onperiod, object=l_fit, original_data=l_x_an_origin, MoreArgs=list("p"=p, "l_time"=l_time), SIMPLIFY=FALSE)
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("x_all", "x_nat", "x_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_allq2 <- mapply(rbinding, l_allq, l_x_an_origin, SIMPLIFY=FALSE)
  a_allq <- array(unlist(l_allq2), dim=c(dim(l_allq2[[1]]), length(l_allq)))
  dimnames(a_allq)[1:2] <- dimnames(l_allq2[[1]])
  dimnames(a_allq)[[3]] <- 1:length(l_allq)
  a_allq
}


imput_aurel_byyear <- function(RR){
  l_na <- which(is.na(RR))
  n_na <- length(l_na)
  print(n_na/length(RR))
  n_inf <- floor(n_na/2)
  n_sup <- n_na - n_inf
  if(n_inf > 0) {
    i_inf <- l_na[1:n_inf]
    RR[i_inf]  <-  0
  }
  if(n_sup > 0) {
    i_sup <- l_na[(n_inf+1) : n_na]
    RR[i_sup]  <-  Inf
  }
  RR
}


# Imputs value in the case of Relative Risk(RR) being equal to NA 
# 
# 
# \code{imput_aure} Imputs value in the case of Relative Risk(RR) being equal
# to NA for each bootstrap sample and for each year.
#
# @param boot_res_RR a matrix of the bootstrap estimate of RR for each year and
# each bootstrap sample
# @return returns a matrix the NA value replaced half of the time by +infinity
# and half by the time by zero
# @examples
# # Needs to be done !
#' @export
imput_aurel <- function(boot_res_RR){
  RR <- apply(boot_res_RR, 1, imput_aurel_byyear)
  t(RR)
}

# Adds a paramater to the array containg the bootstrap results  .
# 
# \code{add_param} Adds a parameter to the the array containing the
# bootstrap results. The added parameters if computed as a function of the
# parameters already presents.

#' Adds a parameter to the the array containing the
#' bootstrap results. The added parameters if computed as a function of the
#' parameters amready presents. The function uses non-standard evaluation to
#' compute the new parameters. For instance, the FAR could be computed from the
#' probabilities p_all and p_nat.
#'
#' @param b_onperiod array containing the estimates from the bootstrap samples.
#' @param operation an expression telling how to compute the new parameter from
#' the others parameters.
#' @param name the name of the new parameters as it appears in the boootstrap
#' array
#' @return returns the bootstrap array with the additional row containg the new
#' parameters computed
#' @examples
#' ans <- compute_far.default(model="cnrm", y="eur_tas", x="gbl_tas", time="year", xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' # get bootstrap samples of p_all p_nat and p_ant
#' bp <- ans$allp
#' # compute the evolution of p_all relatively to its value in 1850
#' bp  <- add_param(bp, operation=p_all/p_all[time==1850], name="p_all_rel")
#' # create the 0.95 confidence intervals data.frame
#' ic  <- get_ic_onperiod(bp, ci=0.95)
#' @export
add_param <- function(b_onperiod, operation, name){
  ans <- array(dim=dim(b_onperiod)+c(1,0,0))
  ans[-dim(ans)[1],,] <- b_onperiod
  rownames_ans <- dimnames(b_onperiod)[[1]]
  dimnames(ans)[2:3] <- dimnames(b_onperiod)[2:3]
  b_onperiod %<>% melt(., varnames=c("param", "time", "bootsample")) %>%
  dcast(bootsample+time~param, data=.)
  b_order  <- order(b_onperiod$time)
  new_param <- eval(substitute(operation), envir=b_onperiod[b_order,])
  rownames_ans %<>% append(., name)
  ans[dim(ans)[1],,][b_order] <- new_param
  dimnames(ans)[[1]] <- rownames_ans
  ans
}

#' Compute confidence intervals from the bootstrap samples
#' 
#' \code{get_ic_onperiod} Compute the confidences intervals for each parameters
#' of the bootstrap samples and at each time.
#' Compute the confidences intervals at confidence level  ci_p for each parameters
#' of the bootstrap samples and at each time. By default, it computes
#' 90%-confidence intervals. The best estimate corresponds to the median of
#' bootstrap estimates.
#' @param b_onperiod array containing the estimates from the bootstrap samples.
#' @param method_name adds a columns method to the results with value given by
#' method_name
#' @param ci_p the level of confidence
#' @param ... additional parameters for the aggregate function use to compute
#' the bootstrap confidence intervales
#' array
#' @return returns a data.frame with columns, method, param, IC_inf, Estim (best
#' estimate), IC_sup .
#' @examples
#' ans <- compute_far.default(model="cnrm", y="eur_tas", x="gbl_tas", time="year", xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' # get bootstrap samples of p_all p_nat and p_ant
#' bp <- ans$allp
#' # compute the evolution of p_all relatively to its value in 1850
#' bp  <- add_param(bp, operation=p_all/p_all[time==1850], name="p_all_rel")
#' # create the 0.95 confidence intervals data.frame
#' ic  <- get_ic_onperiod(bp, ci=0.95)
#' # To be done
#' @export
get_ic_onperiod <- function(b_onperiod, method_name=NULL,  ci_p=0.95, ...){
  alpha <- 1-ci_p
  if (is.null(method_name)) method_name <- deparse(substitute(b_onperiod))
  b_onperiod %<>% melt(., varnames=c("param", "time", "bootsample"))%>% 
  cbind(., "method" = method_name) %>%
  aggregate(value~param+time+method, data=., FUN=quantile, probs=c(alpha/2, 0.5, 1-alpha/2), ...)
  ic <- b_onperiod$value
  colnames(ic) <- c("IC_inf", "Estim", "IC_sup")
  b_onperiod$value <- NULL
  ans <- cbind(b_onperiod, ic)
  attr(ans, "ci_p") = ci_p
  ans
}


