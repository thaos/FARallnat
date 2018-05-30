# Creates function to co compute the FAR 
#
# \code{constrained_far} compute the far for wich the non-stationnary
# components from model fit_o are constrained by the statistical model fit_m.
# i.e. fit_o must have the same fitted parameters has fit_m with the exception
# of the intercept parameters.
# 
# @param far_o result from a function compute_far
# @param far_m result from a function compute_far. It should have the same
# stastical model for y as far_o. the fitted parameters with the exception of
# the intercepts will constrained the statistical model that will be refitted
# for far_o.
# @return the far_o object where the the fitted model far_o$mfit has been
# constrain by the results of far_m
# @examples
# far_cnrm <- compute_far.default("cnrm", y="eur_tas", x="gbl_tas", time="year")
# far_obs <- compute_far.default("obs", y="eur_tas", x="gbl_tas", time="year")
# far_cstr <- constrained_far(far_obs, far_cnrm, model="cnrm")
# @export
constrained_far <- function(far_o, far_m, model, xp=1.6, ci_p=0.9){
  fit_o  <- far_o$mfit
  fit_m  <- far_m$mfit
  fit_c <- mapply(refit, fit_o, fit_m, SIMPLIFY=FALSE)
  dx <- far_o$dx
  allp <- boot_allp_onperiod(fit_c, dx$l_x_an_origin, l_time=sort(unique(dx$l_x_an_origin[[1]]$time)), xp=xp)
  # computing the Relative Risk(RR) for each bootstrap sample
  far <- add_param(allp, operation=p_all/p_nat, name="RR")
  # if RR=NA remplace by either +inf or 0
  far["RR",,] %<>% imput_aurel(.) 
  # add rescaled far (notably, atan(log(RR)) the FAR plot 
  far %<>%  add_param(., operation=al_trans(RR), name="alRR") %>% 
  add_param(., operation=el_trans(RR), name="elRR")
  # compute  confidence intervals for each parameters
  ic_far <- get_ic_onperiod(far, method_name=model, ci_p=ci_p)  
  ic_far$param <- as.character(ic_far$param)
  # ----------------------------------------------------------------------------------
  # return a list with notably all bootsrap decomposition dx and all
  # statistical fits for y.
  far_o$mfit <- fit_c
  far_o$allp <- allp
  far_o$ic_far <- ic_far
  far_o
}
