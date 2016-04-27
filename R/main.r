#' Compute the FAR for all the models
#'
#' \code{compute_forall} computes the FAR evolution for all models according to the compute_far function taken as
#' an argument. if the argument p is provided, it also computes the
#' evolution of the corresponding quantile. 
#'
#' @param compute_far a function used to compute the FAR 
#' @param l_models a character vectors giving the list of models for which the
#' FAR needs to be computed. by default all models present in this packages are
#' used
#' @param xp the selected threshold used to define the FAR
#' @param p the selected probability of exceedance to define the corresponding
#' quantile
#' @param y the name of the variable that will be used as y in the compute_far
#' function
#' @param x the name of the variable that will be used as the covariate x in the compute_far
#' function
#' @param time the name of the variable that will be used as the time variable  in the compute_far
#' function
#' @param res_folder the folder where the results of compute_far for each model
#' are saved as .rds. By default, creates a folder the name of which is the name
#' of the function compute_far
#' @param ... additional parameters needed in the compute_far function if
#' required
#' @param ci_p the level of the confidence intervals
#' @return return a dataframe with the confidence intervals for the FAR and all
#' the other computed parameters (e.g, p and q) at each time
#' @examples
#' # compute the FAR for the CNRM and the IPSL GCMs using a  gam decomposition
#' # and a gaussian fit with only three bootstrap samples
#' ans <- compute_forall(compute_far.default, stat_model=gauss_fit, p=0.01, R=3)
#' # same with ebm decompositon and a gpd fit over the 90% quantil
#' ans <- compute_forall(compute_far.dx_ebm_fit, l_models=c("cnrm", "ipsl"), stat_model=gpd_fit, qthreshold =0.90, p=0.01, R=3)
#' @export
compute_forall <- function(compute_far, l_models=NULL, xp=1.6, p=NULL,  y="eur_tas", x="gbl_tas", time="year", ci_p=0.9, res_folder=deparse(substitute(compute_far)), ...){
  dir.create(res_folder, recursive=TRUE)
  if(is.null(l_models)){
    l_models <- c("bcc", "bnu", "cccma", "cmcc", "cnrm", "csiro", "fio", "gfdl", "giss", "iap", "ichec", "ingv", "inm", "ipsl", "miroc", "mohc", "mpim", "mri", "ncar", "ncc", "obs")
  }
  l_far_files <- paste(res_folder, "/", l_models, ".rds", sep="")
  l_far <- character(length(l_models))
  for(i in seq_along(l_models)){
    cat("------------------------------------\n", l_models[i])
    ans <- compute_far(l_models[i], xp=xp, y=y, x=x, time=time, ci_p=ci_p, ...) 
    if(!is.null(p)) ans <- continue_with_q(ans, l_models[i], p)
    saveRDS(ans, file=l_far_files[i]) 
  }
  datas <- do.call(rbind, mapply(function(model, file){cbind(readRDS(file)$data, method=model)}, model=l_models, file=l_far_files, SIMPLIFY=FALSE))
  ic_fars <- do.call(rbind, mapply(function(model, file){readRDS(file)$ic_far}, model=l_models, file=l_far_files, SIMPLIFY=FALSE))
  merged_res <- merge(ic_fars, datas, by.x=c("time", "method"), by.y=c("year", "method"))
  attr(merged_res, "match.call") <- mget(names(formals()),sys.frame(sys.nframe()))
  merged_res
}

#' Compute the FAR from observations using constraints 
#'
#' \code{compute_forall} Compute the FAR from the observations using the GCM to
#' constrain the evolution of the statistical distribution of the variable of
#' interest y. The relationship between the variable y and the covariate x in the
#' observational should be in the observation the same as in the GCM.
#'
#' @param merged_res a result from the compute_forall function. the dataset
#' "obs" should have been treated in the compute_forall function.
#' @return return a dataframe with the confidence intervals for the constrained FAR and all
#' the other computed parameters (e.g, p and q) at each time.
#' @examples
#' # compute the FAR for the CNRM and the IPSL GCMs using a  gam decomposition
#' # and a gaussian fit with only three bootstrap samples
#' ans <- compute_forall(compute_far.default, stat_model=gauss_fit, p=0.01, R=3)
#' ans_cstr <- continue_with_constrained_far(ans) 
#' summary_plot(ans_cstr)
#' @export
continue_with_constrained_far <- function(merged_res){
  mcall <- attr(merged_res, "match.call")
  x <- mcall$x
  y <- mcall$y
  xp <- mcall$xp
  p <- mcall$p
  ci_p <- mcall$ci_p
  res_folder <- mcall$res_folder
  l_models <- mcall$l_models
  l_models <- l_models[l_models != "obs"]
  l_far_files <- paste(res_folder, "/", l_models, "_cstr.rds", sep="")
  far_o <- readRDS(paste(res_folder,"/obs.rds", sep="")) 
  for(i in seq_along(l_models)){
    cat("------------------------------------\n", l_models[i])
    far_m <- readRDS(paste(res_folder,"/", l_models[i], ".rds", sep="")) 
    ans <- constrained_far(far_o, far_m, model=l_models[i], xp=xp, ci_p=ci_p) 
    if(!is.null(p)) ans <- continue_with_q(ans, l_models[i], p)
    saveRDS(ans, file=l_far_files[i]) 
  }
  datas <- do.call(rbind, mapply(function(model, file){cbind(readRDS(file)$data, method=model)}, model=l_models, file=l_far_files, SIMPLIFY=FALSE))
  ic_fars <- do.call(rbind, mapply(function(model, file){readRDS(file)$ic_far}, model=l_models, file=l_far_files, SIMPLIFY=FALSE))
  merged_res <- merge(ic_fars, datas, by.x=c("time", "method"), by.y=c("year", "method"))
  attr(merged_res, "match.call") <- mcall
  merged_res
}
#' ans <- compute_forall(compute_far.default, stat_model=gauss_fit, p=0.01, R=3)
#' ans_cstr <- continue_with_constrained_far(ans) 
#' summary_plot(ans_cstr)

#' Compute the quantile corresponding to a probability of exceedance p 
#'
#' \code{continue_with_q} computes the quantile corresponding to a probability
#' of exceedance p from the results of a compute_far function
#'
#' @param ans results of a compute_far function
#' @param model the name of the model for which compute_far has been run
#' @param p the probability used to define the quantile
#' @return a list with the same objects as a compute_far function with in
#' addition an object allq which contains all the bootstrap estimates of q
#' @examples
#' # compute the FAR for the CNRM model using a  gam decomposition
#' # and a gaussian fit with only three bootstrap samples
#' ans <- compute_far.default(model="cnrm", y="eur_tas", x="gbl_tas", time="year", xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' ans <- continue_with_q(ans, model="cnrm", p=0.01)
#' @export
continue_with_q <- function(ans, model,  p){
    ci_p <- attr(ans$ic_far, "ci_p")
    allq <- boot_allq_onperiod(ans$mfit, ans$dx$l_x_an_origin, l_time=sort(unique(ans$dx$l_x_an_origin[[1]]$time)), p=p)
    ic_far <- get_ic_onperiod(allq, method_name=model, ci_p=ci_p)  
    ic_far <- unique(rbind(ans$ic_far, ic_far))
    ans$allq <- allq
    ans$ic_far <- ic_far
    ans
}

#' Summary plots of the estimated FAR
#'
#' \code{summary_plot} computes a fex plots  to summarize the results of 
#' compute_forall
#'
#' @param merged_res results of the compute_far_fonction
#' @param pdf_name the name of the pdf file containing the plots. the pdf can be
#' found in the same folder as the results foolders of merged_res
#' @return saves a pdf files called "summary.pdf" with the plots in same folder
#' as the res_forlder of the compute_forall function
#' @examples
#' # compute the FAR for the CNRM and the IPSL GCMs using a  gam decomposition
#' # and a gaussian fit with only three bootstrap samples
#' ans <- compute_forall(compute_far.default, stat_model=gauss_fit, p=0.01, R=3)
#' # summary plot
#' summary_plot(ans)
#' @export
summary_plot <- function(merged_res, pdf_name="summary.pdf"){
  mcall <- attr(merged_res, "match.call")
  x <- mcall$x
  y <- mcall$y
  xp <- mcall$xp
  res_folder <- mcall$res_folder
  col=c(gg_color_hue(length(unique(merged_res$method))-1), "black")
  #######################################################################################################
  default_plot <- function(subset){
    ggplot(data=subset, aes(x=time))+
    # ggtitle(expression(paste("Covariate Estimation : ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
    geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
    geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, fill=param), alpha=0.4) + 
    # facet_grid(method~param) +
    facet_wrap(~method) +
    coord_cartesian(xlim=c(1850,2100))+
    # scale_color_manual(values=col)+
    # scale_fill_manual(values=col)+
    theme(legend.position = "bottom")
  }
  reorder_layer  <- function(gp){
    layers <- gp$layer
    nb_layers <- length(layers)
    if(nb_layers >1) layers <- layers[c(2:nb_layers, 1)]
    gp$layers <- layers
    gp
  }
  p_dx <- default_plot(subset=subset(merged_res, param %in% c("x_all", "x_ant", "x_nat"))) +
  geom_hline(aes_q(yintercept=xp), linetype=2)+
  geom_point(aes_string(y=x), color="black",size=0.9, alpha=0.1) 
  p_dx <- reorder_layer(p_dx)
  p_mu <- default_plot(subset=subset(merged_res, param %in% c("mu_all", "mu_ant", "mu_nat", "threshold_all", "threshold_ant", "threshold_nat"))) +
  geom_hline(aes_q(yintercept=xp), linetype=2)+
  geom_point(aes_string(y=y), color="black",size=0.9, alpha=0.1) 
  p_mu <- reorder_layer(p_mu)
  p_p <- default_plot(subset=subset(merged_res, param %in% c("p_all", "p_ant", "p_nat"))) +
  coord_trans(y="sqrt")
  subset_q <- subset(merged_res, param %in% c("q_all", "q_ant", "q_nat"))
  if(nrow(subset_q) != 0) {
    p_q <- default_plot(subset=subset(merged_res, param %in% c("q_all", "q_ant", "q_nat"))) +
    geom_hline(aes_q(yintercept=xp), linetype=2)+
    geom_point(aes_string(y=y), color="black",size=0.9, alpha=0.1) 
    p_q <- reorder_layer(p_q)
  }
  ########################################################################################################
  p_far <- plot_pannel_far(merged_res, axis_trans="al", main="FAR from 1850 to 2100", col=col)
  pdf(paste(res_folder,"/", pdf_name, sep=""), height=5+round(0.5*length(col)), width=5+round(0.5*length(col)))
      plot(p_dx)
      plot(p_mu)
      plot(p_p)
      plot(p_q)
      plot(p_far)
      plot_far(merged_res, axis_trans="al", main="FAR from 1850 to 2100", col=col)
  dev.off()
}

