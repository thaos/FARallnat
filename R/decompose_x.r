#' Function to decompose the covariate x into an ALL, an ANT and a NAT
#' component 
#'
#' \code{dx.gam_allnat} uses the following  gam model x_all = beta_nat * nat + s(ant)
#' @param bsamples : a list of data.frame with the original dataset merged with the
#' bootstrap EBM simulations
#' @param bindexes : a list of vectors of indexes indicating the rows that will be
#' used for the bootstrap
#' @param y the name of the variable of interest  
#' @param x the name of the variable the covariate that will be decomposed
#' @param all the name of the variable that will be used as the ALL variable
#' @param ant the name of the variable that will be used as the ANT variable
#' @param ghg the name of the variable that will be used as the GHG variable
#' @param others the name of the variable that will be used as the OTHERS variable
#' @param nat the name of the variable that will be used as the NAT variable
#' @param time the name of the variable that will be used as the time
#' @param hnat the name of the variable that indicates wheters the data are from
#' an historical only nat run.
#' @param knots the numbers of knots wanted for the spline in the gam fit
#' @param fixed wether the numbers of knots is allowed to vary or not
#' @param correct_ant_bias whether the variables that are anthropogenic need to
#' have their mean shifted to zero in the pre-industrial.
#' @param pre_ind period of the pre-industrial if the anthropogenic variables need to be
#' shifted.
#' @param reuse_ant_bias wheter only the shift computed on the first bootstrap
#' samples is used to reshift the other bootstrap samples. Otherwise each
#' bootstrap samples is correctied by its own shift.
#' @return two list :
#' \itemize{ 
#' \item bsamples a list of data.frame containing the bootstrap samples with the
#' decomposition of the covariate x into x_all, x_nat and x_ant
#' \item osamples a list of data.frame containing the completed/original
#' datasets for each simulation of the EBM as well as the decomposition of the
#' covariate x into x_all, x_nat and x_ant
#' }
#' @export
dx.gam_allnat <- function(bsamples, bindexes, y="y", x="x", ant="ant", nat="nat", time="time", hnat="hnat", knots=NULL, fixed=FALSE, correct_ant_bias=TRUE, pre_ind=c(1850, 1879), reuse_ant_bias=FALSE){
  osamples <- lapply(bsamples, function(data) data.frame("y" = data[, y], "x" = data[, x], "ant" = data[, ant], "nat" = data[, nat], "time" = data[, time], "hnat" = data[, hnat]))
  psamples <- prepare_osamples(osamples, bindexes)
  if(correct_ant_bias){
    psamples <- boot_correct_ant_bias(psamples, ant_var="ant",  pre_ind=pre_ind, reuse_ant_bias=reuse_ant_bias)
  }
  dx <- boot_gam_allnat(psamples)
}

#' \code{dx.raw} keeps the simulated EBM forcings as it is
#' @rdname dx.gam_allnat 
#' @export
dx.raw <- function(bsamples, bindexes, y="y", all="all", ant="ant", nat="nat", time="time", hnat="hnat"){
  osamples <- lapply(bsamples, function(data) data.frame("y" = data[, y], "x_all" = data[, all], "x_ant" = data[, ant], "x_nat" = data[, nat], "time" = data[, time], "hnat" = data[, hnat]))
  psamples <- prepare_osamples(osamples, bindexes)
  names(psamples) <- c("l_x_an", "l_x_an_origin")
  psamples
}

#' \code{dx.lm_gno} uses the follwing linear model : x_all = beta_nat * nat + beta_ghg * ghg + beta_others * others
#' @export
#' @rdname dx.gam_allnat 
dx.lm_gno <- function(bsamples, bindexes, y="y", x="x", ghg="ghg", nat="nat", others="others", time="time", hnat="hnat", correct_ant_bias=TRUE, pre_ind=c(1850, 1879), reuse_ant_bias=FALSE){
  osamples <- lapply(bsamples, function(data) data.frame("y" = data[, y], "x" = data[, x], "ghg" = data[, ghg], "nat" = data[, nat], "others"=data[, others], "time" = data[, time], "hnat" = data[, hnat]))
  psamples <- prepare_osamples(osamples, bindexes)
  if(correct_ant_bias){
    psamples <- boot_correct_ant_bias(psamples, ant_var=c("ghg", "others"),  pre_ind=pre_ind, reuse_ant_bias=reuse_ant_bias)
  }
  dx <- boot_lm_gno(psamples)
}

#' \code{dx.gam_gno} uses the follwing linear model : gam model x_all = beta_nat * nat + beta_ghg * ghg + s(others)
#' @export
#' @rdname dx.gam_allnat 
dx.gam_gno <- function(bsamples, bindexes, y="y", x="x", ghg="ghg", nat="nat", others="others", time="time", hnat="hnat", correct_ant_bias=TRUE, pre_ind=c(1850, 1879), reuse_ant_bias=FALSE){
  osamples <- lapply(bsamples, function(data) data.frame("y" = data[, y], "x" = data[, x], "ghg" = data[, ghg], "nat" = data[, nat], "others"=data[, others], "time" = data[, time], "hnat" = data[, hnat]))
  psamples <- prepare_osamples(osamples, bindexes)
  if(correct_ant_bias){
    psamples <- boot_correct_ant_bias(psamples, ant_var=c("ghg", "others"),  pre_ind=pre_ind, reuse_ant_bias=reuse_ant_bias)
  }
  dx <- boot_gam_gno(psamples)
}

#' \code{dx.gam_gno} fits EBM using ordinary least squares for the decomposition
#' @export
#' @rdname dx.gam_allnat 
dx.ebm_fit <- function(bsamples, bindexes, y="y", x="x", time="time", hnat="hnat"){
  osamples <- lapply(bsamples, function(data) data.frame("y" = data[, y], "x" = data[, x], "time" = data[, time], "hnat" = data[, hnat]))
  psamples <- prepare_osamples(osamples, bindexes)
  dx <- boot_ebm_allnat(psamples)
}
